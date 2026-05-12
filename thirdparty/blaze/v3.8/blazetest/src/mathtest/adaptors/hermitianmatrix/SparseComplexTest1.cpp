//=================================================================================================
/*!
//  \file src/mathtest/adaptors/hermitianmatrix/SparseComplexTest1.cpp
//  \brief Source file for the HermitianMatrix sparse complex test (part 1)
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
#include <blaze/math/StaticMatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/adaptors/hermitianmatrix/SparseComplexTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace hermitianmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the HermitianMatrix sparse complex test.
//
// \exception std::runtime_error Operation error detected.
*/
SparseComplexTest::SparseComplexTest()
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testSchurAssign();
   testMultAssign();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the HermitianMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the HermitianMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   // Default constructor (CompressedMatrix)
   {
      test_ = "Row-major HermitianMatrix default constructor (CompressedMatrix)";

      const HT herm;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   // Size constructor (CompressedMatrix)
   {
      test_ = "Row-major HermitianMatrix size constructor (CompressedMatrix)";

      const HT herm( 2UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkNonZeros( herm, 0UL );
   }


   //=====================================================================================
   // Row-major list initialization
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Row-major HermitianMatrix initializer list constructor (complete list)";

      const HT herm{ { cplx(1, 0), cplx(2, 2), cplx(3,-3) },
                     { cplx(2,-2), cplx(4, 0), cplx(0, 5) },
                     { cplx(3, 3), cplx(0,-5), cplx(6, 0) } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 9UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 2) || herm(0,2) != cplx(3,-3) ||
          herm(1,0) != cplx(2,-2) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 5) ||
          herm(2,0) != cplx(3, 3) || herm(2,1) != cplx(0,-5) || herm(2,2) != cplx(6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 2) (3,-3) )\n"
                                     "( (2,-2) (4, 0) (0, 5) )\n"
                                     "( (3, 3) (0,-5) (6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Row-major HermitianMatrix initializer list constructor (incomplete list)";

      const HT herm{ { cplx(1, 0), cplx(2, 2), cplx(3,-3) },
                     { cplx(2,-2), cplx(4, 0) },
                     { cplx(3, 3), cplx(0, 0), cplx(6, 0) } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 2) || herm(0,2) != cplx(3,-3) ||
          herm(1,0) != cplx(2,-2) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 0) ||
          herm(2,0) != cplx(3, 3) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 2) (3,-3) )\n"
                                     "( (2,-2) (4, 0) (0, 0) )\n"
                                     "( (3, 3) (0, 0) (6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Row-major HermitianMatrix copy constructor (0x0)";

      const HT herm1;
      const HT herm2( herm1 );

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Row-major HermitianMatrix copy constructor (3x3)";

      HT herm1( 3UL );
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      const HT herm2( herm1 );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Row-major HermitianMatrix move constructor (0x0)";

      HT herm1;
      HT herm2( std::move( herm1 ) );

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Row-major HermitianMatrix move constructor (3x3)";

      HT herm1( 3UL );
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      HT herm2( std::move( herm1 ) );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Row-major HermitianMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<cplx,blaze::rowMajor> mat;
      const HT herm( mat );

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Conversion constructor (Hermitian)
   {
      test_ = "Row-major HermitianMatrix conversion constructor (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      const HT herm( mat );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (invalid diagonal)
   {
      test_ = "Row-major HermitianMatrix conversion constructor (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 1);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      try {
         const HT herm( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-Hermitian HermitianMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (invalid pair)
   {
      test_ = "Row-major HermitianMatrix conversion constructor (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7, 3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      try {
         const HT herm( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-Hermitian HermitianMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (HermitianMatrix)
   {
      test_ = "Row-major HermitianMatrix conversion constructor (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<cplx,3UL,3UL,blaze::rowMajor> > herm1;
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      const HT herm2( herm1 );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   // Default constructor (CompressedMatrix)
   {
      test_ = "Column-major HermitianMatrix default constructor (CompressedMatrix)";

      const OHT herm;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   // Size constructor (CompressedMatrix)
   {
      test_ = "Column-major HermitianMatrix size constructor (CompressedMatrix)";

      const OHT herm( 2UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkNonZeros( herm, 0UL );
   }


   //=====================================================================================
   // Column-major list initialization
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Column-major HermitianMatrix initializer list constructor (complete list)";

      const OHT herm{ { cplx(1, 0), cplx(2, 2), cplx(3,-3) },
                      { cplx(2,-2), cplx(4, 0), cplx(0, 5) },
                      { cplx(3, 3), cplx(0,-5), cplx(6, 0) } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 9UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 2) || herm(0,2) != cplx(3,-3) ||
          herm(1,0) != cplx(2,-2) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 5) ||
          herm(2,0) != cplx(3, 3) || herm(2,1) != cplx(0,-5) || herm(2,2) != cplx(6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 2) (3,-3) )\n"
                                     "( (2,-2) (4, 0) (0, 5) )\n"
                                     "( (3, 3) (0,-5) (6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Column-major HermitianMatrix initializer list constructor (incomplete list)";

      const OHT herm{ { cplx(1, 0), cplx(2, 2), cplx(3,-3) },
                      { cplx(2,-2), cplx(4, 0) },
                      { cplx(3, 3), cplx(0, 0), cplx(6, 0) } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 2) || herm(0,2) != cplx(3,-3) ||
          herm(1,0) != cplx(2,-2) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 0) ||
          herm(2,0) != cplx(3, 3) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 2) (3,-3) )\n"
                                     "( (2,-2) (4, 0) (0, 0) )\n"
                                     "( (3, 3) (0, 0) (6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Column-major HermitianMatrix copy constructor (0x0)";

      const OHT herm1;
      const OHT herm2( herm1 );

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Column-major HermitianMatrix copy constructor (3x3)";

      OHT herm1( 3UL );
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      const OHT herm2( herm1 );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Column-major HermitianMatrix move constructor (0x0)";

      OHT herm1;
      OHT herm2( std::move( herm1 ) );

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Column-major HermitianMatrix move constructor (3x3)";

      OHT herm1( 3UL );
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      OHT herm2( std::move( herm1 ) );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Column-major HermitianMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<cplx,blaze::columnMajor> mat;
      const OHT herm( mat );

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Conversion constructor (Hermitian)
   {
      test_ = "Column-major HermitianMatrix conversion constructor (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      const OHT herm( mat );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (invalid diagonal)
   {
      test_ = "Column-major HermitianMatrix conversion constructor (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 1);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      try {
         const OHT herm( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-Hermitian HermitianMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (invalid pair)
   {
      test_ = "Column-major HermitianMatrix conversion constructor (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7, 3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      try {
         const OHT herm( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-Hermitian HermitianMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (HermitianMatrix)
   {
      test_ = "Column-major HermitianMatrix conversion constructor (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<cplx,3UL,3UL,blaze::columnMajor> > herm1;
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      const OHT herm2( herm1 );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the HermitianMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the HermitianMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testAssignment()
{
   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Row-major HermitianMatrix initializer list assignment (complete list)";

      HT herm;
      herm = { { cplx(1, 0), cplx(2, 2), cplx(3,-3) },
               { cplx(2,-2), cplx(4, 0), cplx(0, 5) },
               { cplx(3, 3), cplx(0,-5), cplx(6, 0) } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 9UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 2) || herm(0,2) != cplx(3,-3) ||
          herm(1,0) != cplx(2,-2) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 5) ||
          herm(2,0) != cplx(3, 3) || herm(2,1) != cplx(0,-5) || herm(2,2) != cplx(6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 2) (3,-3) )\n"
                                     "( (2,-2) (4, 0) (0, 5) )\n"
                                     "( (3, 3) (0,-5) (6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Row-major HermitianMatrix initializer list assignment (incomplete list)";

      HT herm;
      herm = { { cplx(1, 0), cplx(2, 2), cplx(3,-3) },
               { cplx(2,-2), cplx(4, 0) },
               { cplx(3, 3), cplx(0, 0), cplx(6, 0) } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 2) || herm(0,2) != cplx(3,-3) ||
          herm(1,0) != cplx(2,-2) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 0) ||
          herm(2,0) != cplx(3, 3) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 2) (3,-3) )\n"
                                     "( (2,-2) (4, 0) (0, 0) )\n"
                                     "( (3, 3) (0, 0) (6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Row-major HermitianMatrix copy assignment (0x0)";

      HT herm1, herm2;

      herm2 = herm1;

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Row-major HermitianMatrix copy assignment (3x3)";

      HT herm1( 3UL );
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      HT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Row-major HermitianMatrix move assignment (0x0)";

      HT herm1, herm2;

      herm2 = std::move( herm1 );

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Row-major HermitianMatrix move assignment (3x3)";

      HT herm1( 3UL );
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      HT herm2;
      herm2 = std::move( herm1 );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major HermitianMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<cplx,blaze::rowMajor> mat;

      HT herm;
      herm = mat;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Row-major/row-major dense matrix assignment (Hermitian)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      HT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (Hermitian)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      HT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix assignment (invalid diagonal)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 1);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      try {
         HT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix assignment (invalid diagonal)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 1);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      try {
         HT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix assignment (invalid pair)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7, 3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      try {
         HT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix assignment (invalid pair)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7, 3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      try {
         HT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<cplx,3UL,3UL,blaze::rowMajor> > herm1;
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      HT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<cplx,3UL,3UL,blaze::columnMajor> > herm1;
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      HT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major HermitianMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<cplx,blaze::rowMajor> mat;

      HT herm;
      herm = mat;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Row-major/row-major sparse matrix assignment (Hermitian)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 8UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,2) = cplx( 3, 0);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      HT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 8UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (Hermitian)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 8UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,2) = cplx( 3, 0);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      HT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 8UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix assignment (invalid diagonal)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 1);
      mat(2,0) = cplx( 7,-3);
      mat(2,2) = cplx( 3, 0);

      try {
         HT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix assignment (invalid diagonal)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 1);
      mat(2,0) = cplx( 7,-3);
      mat(2,2) = cplx( 3, 0);

      try {
         HT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix assignment (invalid pair)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(2,0) = cplx( 7, 3);
      mat(2,2) = cplx( 3, 0);

      try {
         HT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix assignment (invalid pair)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(2,0) = cplx( 7, 3);
      mat(2,2) = cplx( 3, 0);

      try {
         HT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<blaze::complex<long>,blaze::rowMajor> > herm1( 3UL, 7UL );
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      HT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<cplx,blaze::columnMajor> > herm1( 3UL, 7UL );
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      HT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Column-major HermitianMatrix initializer list assignment (complete list)";

      OHT herm;
      herm = { { cplx(1, 0), cplx(2, 2), cplx(3,-3) },
               { cplx(2,-2), cplx(4, 0), cplx(0, 5) },
               { cplx(3, 3), cplx(0,-5), cplx(6, 0) } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 9UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 2) || herm(0,2) != cplx(3,-3) ||
          herm(1,0) != cplx(2,-2) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 5) ||
          herm(2,0) != cplx(3, 3) || herm(2,1) != cplx(0,-5) || herm(2,2) != cplx(6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 2) (3,-3) )\n"
                                     "( (2,-2) (4, 0) (0, 5) )\n"
                                     "( (3, 3) (0,-5) (6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Column-major HermitianMatrix initializer list assignment (incomplete list)";

      OHT herm;
      herm = { { cplx(1, 0), cplx(2, 2), cplx(3,-3) },
               { cplx(2,-2), cplx(4, 0) },
               { cplx(3, 3), cplx(0, 0), cplx(6, 0) } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 2) || herm(0,2) != cplx(3,-3) ||
          herm(1,0) != cplx(2,-2) || herm(1,1) != cplx(4, 0) || herm(1,2) != cplx(0, 0) ||
          herm(2,0) != cplx(3, 3) || herm(2,1) != cplx(0, 0) || herm(2,2) != cplx(6, 0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 2) (3,-3) )\n"
                                     "( (2,-2) (4, 0) (0, 0) )\n"
                                     "( (3, 3) (0, 0) (6, 0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Column-major HermitianMatrix copy assignment (0x0)";

      OHT herm1, herm2;

      herm2 = herm1;

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Column-major HermitianMatrix copy assignment (3x3)";

      OHT herm1( 3UL );
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      OHT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Column-major HermitianMatrix move assignment (0x0)";

      OHT herm1, herm2;

      herm2 = std::move( herm1 );

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Column-major HermitianMatrix move assignment (3x3)";

      OHT herm1( 3UL );
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      OHT herm2;
      herm2 = std::move( herm1 );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major HermitianMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<cplx,blaze::rowMajor> mat;

      OHT herm;
      herm = mat;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Column-major/row-major dense matrix assignment (Hermitian)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      OHT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (Hermitian)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      OHT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix assignment (invalid diagonal)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 1);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      try {
         OHT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix assignment (invalid diagonal)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 1);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      try {
         OHT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix assignment (invalid pair)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7, 3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      try {
         OHT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix assignment (invalid pair)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(1,2) = cplx( 0, 0);
      mat(2,0) = cplx( 7, 3);
      mat(2,1) = cplx( 0, 0);
      mat(2,2) = cplx( 3, 0);

      try {
         OHT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<cplx,3UL,3UL,blaze::rowMajor> > herm1;
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      OHT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<cplx,3UL,3UL,blaze::columnMajor> > herm1;
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      OHT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major HermitianMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<cplx,blaze::rowMajor> mat;

      OHT herm;
      herm = mat;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Column-major/row-major sparse matrix assignment (Hermitian)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 8UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,2) = cplx( 3, 0);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      OHT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 8UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (Hermitian)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 8UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(2,0) = cplx( 7,-3);
      mat(2,2) = cplx( 3, 0);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      OHT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 8UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix assignment (invalid diagonal)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 1);
      mat(2,0) = cplx( 7,-3);
      mat(2,2) = cplx( 3, 0);

      try {
         OHT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix assignment (invalid diagonal)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 1);
      mat(2,0) = cplx( 7,-3);
      mat(2,2) = cplx( 3, 0);

      try {
         OHT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix assignment (invalid pair)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(2,0) = cplx( 7, 3);
      mat(2,2) = cplx( 3, 0);

      try {
         OHT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix assignment (invalid pair)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = cplx( 1, 0);
      mat(0,1) = cplx(-4,-1);
      mat(0,2) = cplx( 7, 3);
      mat(1,0) = cplx(-4, 1);
      mat(1,1) = cplx( 2, 0);
      mat(2,0) = cplx( 7, 3);
      mat(2,2) = cplx( 3, 0);

      try {
         OHT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<cplx,blaze::rowMajor> > herm1( 3UL, 7UL );
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      OHT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<blaze::complex<long>,blaze::columnMajor> > herm1( 3UL, 7UL );
      herm1(0,0) = cplx( 1, 0);
      herm1(0,1) = cplx(-4,-1);
      herm1(0,2) = cplx( 7, 3);
      herm1(1,1) = cplx( 2, 0);
      herm1(2,2) = cplx( 3, 0);

      OHT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-4,-1) || herm2(0,2) != cplx(7,3) ||
          herm2(1,0) != cplx(-4, 1) || herm2(1,1) != cplx( 2, 0) || herm2(1,2) != cplx(0,0) ||
          herm2(2,0) != cplx( 7,-3) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-4,-1) (7,3) )\n"
                                     "( (-4, 1) ( 2, 0) (0,0) )\n"
                                     "( ( 7,-3) ( 0, 0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the HermitianMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   // Row-major/row-major dense matrix addition assignment (Hermitian)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix addition assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (Hermitian)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix addition assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix addition assignment (invalid diagonal)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix addition assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix addition assignment (invalid diagonal)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix addition assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix addition assignment (invalid pair)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix addition assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix addition assignment (invalid pair)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix addition assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix addition assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix addition assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::rowMajor> > herm1( 3UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix addition assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::columnMajor> > herm1( 3UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix addition assignment (Hermitian)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix addition assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 8UL );
      checkNonZeros( herm, 8UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 3UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix addition assignment (Hermitian)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix addition assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 8UL );
      checkNonZeros( herm, 8UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 3UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix addition assignment (invalid diagonal)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix addition assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix addition assignment (invalid diagonal)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix addition assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix addition assignment (invalid pair)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix addition assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix addition assignment (invalid pair)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix addition assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix addition assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix addition assignment (HermitianMatrix)";

      HT herm1( 3UL, 5UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix addition assignment (HermitianMatrix)";

      OHT herm1( 3UL, 5UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   // Column-major/row-major dense matrix addition assignment (Hermitian)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix addition assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (Hermitian)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix addition assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix addition assignment (invalid diagonal)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix addition assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix addition assignment (invalid diagonal)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix addition assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix addition assignment (invalid pair)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix addition assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix addition assignment (invalid pair)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix addition assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix addition assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix addition assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::rowMajor> > herm1( 3UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix addition assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::columnMajor> > herm1( 3UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix addition assignment (Hermitian)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix addition assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 8UL );
      checkNonZeros( herm, 8UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix addition assignment (Hermitian)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix addition assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 8UL );
      checkNonZeros( herm, 8UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix addition assignment (invalid diagonal)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix addition assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix addition assignment (invalid diagonal)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix addition assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix addition assignment (invalid pair)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix addition assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix addition assignment (invalid pair)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix addition assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix addition assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix addition assignment (HermitianMatrix)";

      HT herm1( 3UL, 5UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix addition assignment (HermitianMatrix)";

      OHT herm1( 3UL, 5UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the HermitianMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major dense matrix subtraction assignment (Hermitian)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix subtraction assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (Hermitian)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix subtraction assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix subtraction assignment (invalid diagonal)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix subtraction assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3,-9);
      mat(2,0) = cplx(-6,-2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix subtraction assignment (invalid diagonal)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix subtraction assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3,-9);
      mat(2,0) = cplx(-6,-2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix subtraction assignment (invalid pair)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix subtraction assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 8);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-9);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix subtraction assignment (invalid pair)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix subtraction assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 8);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-9);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix subtraction assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::rowMajor> > herm1( 3UL );
      herm1(0,1) = cplx( 2,1);
      herm1(0,2) = cplx(-6,2);
      herm1(1,1) = cplx(-3,0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix subtraction assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::columnMajor> > herm1( 3UL );
      herm1(0,1) = cplx( 2,1);
      herm1(0,2) = cplx(-6,2);
      herm1(1,1) = cplx(-3,0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix subtraction assignment (Hermitian)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix subtraction assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-2);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 8UL );
      checkNonZeros( herm, 8UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 3UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix subtraction assignment (Hermitian)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix subtraction assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-2);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 8UL );
      checkNonZeros( herm, 8UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 3UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix subtraction assignment (invalid diagonal)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix subtraction assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3,-9);
      mat(2,0) = cplx(-6,-2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix subtraction assignment (invalid diagonal)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix subtraction assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3,-9);
      mat(2,0) = cplx(-6,-2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix subtraction assignment (invalid pair)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix subtraction assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 8);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-9);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix subtraction assignment (invalid pair)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix subtraction assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 8);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-9);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix subtraction assignment (HermitianMatrix)";

      HT herm1( 3UL, 5UL );
      herm1(0,1) = cplx( 2,1);
      herm1(0,2) = cplx(-6,2);
      herm1(1,1) = cplx(-3,0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix subtraction assignment (HermitianMatrix)";

      OHT herm1( 3UL, 5UL );
      herm1(0,1) = cplx( 2,1);
      herm1(0,2) = cplx(-6,2);
      herm1(1,1) = cplx(-3,0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major dense matrix subtraction assignment (Hermitian)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix subtraction assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (Hermitian)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix subtraction assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix subtraction assignment (invalid diagonal)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix subtraction assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3,-9);
      mat(2,0) = cplx(-6,-2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix subtraction assignment (invalid diagonal)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix subtraction assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3,-9);
      mat(2,0) = cplx(-6,-2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix subtraction assignment (invalid pair)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix subtraction assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 8);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-9);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix subtraction assignment (invalid pair)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix subtraction assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 8);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-9);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix subtraction assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::rowMajor> > herm1( 3UL );
      herm1(0,1) = cplx( 2,1);
      herm1(0,2) = cplx(-6,2);
      herm1(1,1) = cplx(-3,0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix subtraction assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::columnMajor> > herm1( 3UL );
      herm1(0,1) = cplx( 2,1);
      herm1(0,2) = cplx(-6,2);
      herm1(1,1) = cplx(-3,0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix subtraction assignment (Hermitian)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix subtraction assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-2);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 8UL );
      checkNonZeros( herm, 8UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix subtraction assignment (Hermitian)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix subtraction assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-2);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 8UL );
      checkNonZeros( herm, 8UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 3UL );

      if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-6,-2) || herm(0,2) != cplx(13,1) ||
          herm(1,0) != cplx(-6, 2) || herm(1,1) != cplx( 5, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(13,-1) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix subtraction assignment (invalid diagonal)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix subtraction assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3,-9);
      mat(2,0) = cplx(-6,-2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix subtraction assignment (invalid diagonal)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix subtraction assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 2);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3,-9);
      mat(2,0) = cplx(-6,-2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix subtraction assignment (invalid pair)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix subtraction assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 8);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-9);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix subtraction assignment (invalid pair)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix subtraction assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx( 2, 1);
      mat(0,2) = cplx(-6, 8);
      mat(1,0) = cplx( 2,-1);
      mat(1,1) = cplx(-3, 0);
      mat(2,0) = cplx(-6,-9);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix subtraction assignment (HermitianMatrix)";

      HT herm1( 3UL, 5UL );
      herm1(0,1) = cplx( 2,1);
      herm1(0,2) = cplx(-6,2);
      herm1(1,1) = cplx(-3,0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix subtraction assignment (HermitianMatrix)";

      OHT herm1( 3UL, 5UL );
      herm1(0,1) = cplx( 2,1);
      herm1(0,2) = cplx(-6,2);
      herm1(1,1) = cplx(-3,0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 1, 0) || herm2(0,1) != cplx(-6,-2) || herm2(0,2) != cplx(13,1) ||
          herm2(1,0) != cplx(-6, 2) || herm2(1,1) != cplx( 5, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(13,-1) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 1, 0) (-6,-2) (13,1) )\n"
                                     "( (-6, 2) ( 5, 0) ( 0,0) )\n"
                                     "( (13,-1) ( 0, 0) ( 3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the HermitianMatrix Schur product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment operators of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testSchurAssign()
{
   //=====================================================================================
   // Row-major dense matrix Schur product assignment
   //=====================================================================================

   // Row-major/row-major dense matrix Schur product assignment (Hermitian)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix Schur product assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 0, 0) || herm(0,1) != cplx(7,6) || herm(0,2) != cplx(48,4) ||
          herm(1,0) != cplx( 7,-6) || herm(1,1) != cplx(6,0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(48,-4) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix Schur product assignment (Hermitian)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix Schur product assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 0, 0) || herm(0,1) != cplx(7,6) || herm(0,2) != cplx(48,4) ||
          herm(1,0) != cplx( 7,-6) || herm(1,1) != cplx(6,0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(48,-4) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix Schur product assignment (invalid diagonal)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix Schur product assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix Schur product assignment (invalid diagonal)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix Schur product assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix Schur product assignment (invalid pair)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix Schur product assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix Schur product assignment (invalid pair)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix Schur product assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix Schur product assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::rowMajor> > herm1( 3UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 0, 0) || herm2(0,1) != cplx(7,6) || herm2(0,2) != cplx(48,4) ||
          herm2(1,0) != cplx( 7,-6) || herm2(1,1) != cplx(6,0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(48,-4) || herm2(2,1) != cplx(0,0) || herm2(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix Schur product assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::columnMajor> > herm1( 3UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 0, 0) || herm2(0,1) != cplx(7,6) || herm2(0,2) != cplx(48,4) ||
          herm2(1,0) != cplx( 7,-6) || herm2(1,1) != cplx(6,0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(48,-4) || herm2(2,1) != cplx(0,0) || herm2(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix Schur product assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix Schur product assignment (Hermitian)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix Schur product assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );

      if( herm(0,0) != cplx( 0, 0) || herm(0,1) != cplx(7,6) || herm(0,2) != cplx(48,4) ||
          herm(1,0) != cplx( 7,-6) || herm(1,1) != cplx(6,0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(48,-4) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix Schur product assignment (Hermitian)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix Schur product assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );

      if( herm(0,0) != cplx( 0, 0) || herm(0,1) != cplx(7,6) || herm(0,2) != cplx(48,4) ||
          herm(1,0) != cplx( 7,-6) || herm(1,1) != cplx(6,0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(48,-4) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix Schur product assignment (invalid diagonal)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix Schur product assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix Schur product assignment (invalid diagonal)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix Schur product assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix Schur product assignment (invalid pair)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix Schur product assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix Schur product assignment (invalid pair)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix Schur product assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix Schur product assignment (HermitianMatrix)";

      HT herm1( 3UL, 5UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 5UL );
      checkNonZeros( herm2, 5UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 1UL );

      if( herm2(0,0) != cplx( 0, 0) || herm2(0,1) != cplx(7,6) || herm2(0,2) != cplx(48,4) ||
          herm2(1,0) != cplx( 7,-6) || herm2(1,1) != cplx(6,0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(48,-4) || herm2(2,1) != cplx(0,0) || herm2(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix Schur product assignment (HermitianMatrix)";

      OHT herm1( 3UL, 5UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 5UL );
      checkNonZeros( herm2, 5UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 1UL );

      if( herm2(0,0) != cplx( 0, 0) || herm2(0,1) != cplx(7,6) || herm2(0,2) != cplx(48,4) ||
          herm2(1,0) != cplx( 7,-6) || herm2(1,1) != cplx(6,0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(48,-4) || herm2(2,1) != cplx(0,0) || herm2(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix Schur product assignment
   //=====================================================================================

   // Column-major/row-major dense matrix Schur product assignment (Hermitian)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix Schur product assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 0, 0) || herm(0,1) != cplx(7,6) || herm(0,2) != cplx(48,4) ||
          herm(1,0) != cplx( 7,-6) || herm(1,1) != cplx(6,0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(48,-4) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix Schur product assignment (Hermitian)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix Schur product assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 0, 0) || herm(0,1) != cplx(7,6) || herm(0,2) != cplx(48,4) ||
          herm(1,0) != cplx( 7,-6) || herm(1,1) != cplx(6,0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(48,-4) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix Schur product assignment (invalid diagonal)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix Schur product assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix Schur product assignment (invalid diagonal)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix Schur product assignment (invalid diagonal)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix Schur product assignment (invalid pair)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix Schur product assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix Schur product assignment (invalid pair)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix Schur product assignment (invalid pair)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix Schur product assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::rowMajor> > herm1( 3UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 0, 0) || herm2(0,1) != cplx(7,6) || herm2(0,2) != cplx(48,4) ||
          herm2(1,0) != cplx( 7,-6) || herm2(1,1) != cplx(6,0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(48,-4) || herm2(2,1) != cplx(0,0) || herm2(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix Schur product assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::columnMajor> > herm1( 3UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 0, 0) || herm2(0,1) != cplx(7,6) || herm2(0,2) != cplx(48,4) ||
          herm2(1,0) != cplx( 7,-6) || herm2(1,1) != cplx(6,0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(48,-4) || herm2(2,1) != cplx(0,0) || herm2(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix Schur product assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix Schur product assignment (Hermitian)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix Schur product assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );

      if( herm(0,0) != cplx( 0, 0) || herm(0,1) != cplx(7,6) || herm(0,2) != cplx(48,4) ||
          herm(1,0) != cplx( 7,-6) || herm(1,1) != cplx(6,0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(48,-4) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix Schur product assignment (Hermitian)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix Schur product assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 2);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 5UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );

      if( herm(0,0) != cplx( 0, 0) || herm(0,1) != cplx(7,6) || herm(0,2) != cplx(48,4) ||
          herm(1,0) != cplx( 7,-6) || herm(1,1) != cplx(6,0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(48,-4) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix Schur product assignment (invalid diagonal)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix Schur product assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix Schur product assignment (invalid diagonal)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix Schur product assignment (invalid diagonal)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-2);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 9);
      mat(2,0) = cplx( 6, 2);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix Schur product assignment (invalid pair)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix Schur product assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix Schur product assignment (invalid pair)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix Schur product assignment (invalid pair)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 5UL );
      mat(0,1) = cplx(-2,-1);
      mat(0,2) = cplx( 6,-8);
      mat(1,0) = cplx(-2, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 9);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix Schur product assignment (HermitianMatrix)";

      HT herm1( 3UL, 5UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 5UL );
      checkNonZeros( herm2, 5UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 1UL );

      if( herm2(0,0) != cplx( 0, 0) || herm2(0,1) != cplx(7,6) || herm2(0,2) != cplx(48,4) ||
          herm2(1,0) != cplx( 7,-6) || herm2(1,1) != cplx(6,0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(48,-4) || herm2(2,1) != cplx(0,0) || herm2(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix Schur product assignment (HermitianMatrix)";

      OHT herm1( 3UL, 5UL );
      herm1(0,1) = cplx(-2,-1);
      herm1(0,2) = cplx( 6,-2);
      herm1(1,1) = cplx( 3, 0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 5UL );
      checkNonZeros( herm2, 5UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 1UL );

      if( herm2(0,0) != cplx( 0, 0) || herm2(0,1) != cplx(7,6) || herm2(0,2) != cplx(48,4) ||
          herm2(1,0) != cplx( 7,-6) || herm2(1,1) != cplx(6,0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(48,-4) || herm2(2,1) != cplx(0,0) || herm2(2,2) != cplx( 0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n( ( 0, 0) (7,6) (48,4) )\n"
                                     "( ( 7,-6) (6,0) ( 0,0) )\n"
                                     "( (48,-4) (0,0) ( 0,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the HermitianMatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseComplexTest::testMultAssign()
{
   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major dense matrix multiplication assignment (Hermitian)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix multiplication assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,0) = cplx(2,0);
      mat(1,1) = cplx(2,0);
      mat(2,2) = cplx(2,0);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm *= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 2, 0) || herm(0,1) != cplx(-8,-2) || herm(0,2) != cplx(14,6) ||
          herm(1,0) != cplx(-8, 2) || herm(1,1) != cplx( 4, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(14,-6) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (Hermitian)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix multiplication assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,0) = cplx(2,0);
      mat(1,1) = cplx(2,0);
      mat(2,2) = cplx(2,0);

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm *= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 2, 0) || herm(0,1) != cplx(-8,-2) || herm(0,2) != cplx(14,6) ||
          herm(1,0) != cplx(-8, 2) || herm(1,1) != cplx( 4, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(14,-6) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix multiplication assignment (non-Hermitian)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix multiplication assignment (non-Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-2);
      mat(0,2) = cplx( 6, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 3);

      HT herm( 3UL );
      herm(0,0) = cplx( 1,0);
      herm(0,1) = cplx(-4,0);
      herm(0,2) = cplx( 7,0);
      herm(1,1) = cplx( 2,0);
      herm(2,2) = cplx( 3,0);

      try {
         herm *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix multiplication assignment (non-Hermitian)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix multiplication assignment (non-Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-2);
      mat(0,2) = cplx( 6, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 3);

      HT herm( 3UL );
      herm(0,0) = cplx( 1,0);
      herm(0,1) = cplx(-4,0);
      herm(0,2) = cplx( 7,0);
      herm(1,1) = cplx( 2,0);
      herm(2,2) = cplx( 3,0);

      try {
         herm *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix multiplication assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix multiplication assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::rowMajor> > herm1( 3UL );
      herm1(0,0) = cplx(2,0);
      herm1(1,1) = cplx(2,0);
      herm1(2,2) = cplx(2,0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 *= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 2, 0) || herm2(0,1) != cplx(-8,-2) || herm2(0,2) != cplx(14,6) ||
          herm2(1,0) != cplx(-8, 2) || herm2(1,1) != cplx( 4, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(14,-6) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix multiplication assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::columnMajor> > herm1( 3UL );
      herm1(0,0) = cplx(2,0);
      herm1(1,1) = cplx(2,0);
      herm1(2,2) = cplx(2,0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 *= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 2, 0) || herm2(0,1) != cplx(-8,-2) || herm2(0,2) != cplx(14,6) ||
          herm2(1,0) != cplx(-8, 2) || herm2(1,1) != cplx( 4, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(14,-6) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix multiplication assignment (Hermitian)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix multiplication assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = cplx(2,0);
      mat(1,1) = cplx(2,0);
      mat(2,2) = cplx(2,0);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm *= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 2, 0) || herm(0,1) != cplx(-8,-2) || herm(0,2) != cplx(14,6) ||
          herm(1,0) != cplx(-8, 2) || herm(1,1) != cplx( 4, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(14,-6) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (Hermitian)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix multiplication assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = cplx(2,0);
      mat(1,1) = cplx(2,0);
      mat(2,2) = cplx(2,0);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      HT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm *= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 2, 0) || herm(0,1) != cplx(-8,-2) || herm(0,2) != cplx(14,6) ||
          herm(1,0) != cplx(-8, 2) || herm(1,1) != cplx( 4, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(14,-6) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix multiplication assignment (non-Hermitian)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix multiplication assignment (non-Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = cplx(-2,-2);
      mat(0,2) = cplx( 6, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 3);

      HT herm( 3UL );
      herm(0,0) = cplx( 1,0);
      herm(0,1) = cplx(-4,0);
      herm(0,2) = cplx( 7,0);
      herm(1,1) = cplx( 2,0);
      herm(2,2) = cplx( 3,0);

      try {
         herm *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix multiplication assignment (non-Hermitian)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix multiplication assignment (non-Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = cplx(-2,-2);
      mat(0,2) = cplx( 6, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 3);

      HT herm( 3UL );
      herm(0,0) = cplx( 1,0);
      herm(0,1) = cplx(-4,0);
      herm(0,2) = cplx( 7,0);
      herm(1,1) = cplx( 2,0);
      herm(2,2) = cplx( 3,0);

      try {
         herm *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix multiplication assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix multiplication assignment (HermitianMatrix)";

      HT herm1( 3UL, 3UL );
      herm1(0,0) = cplx(2,0);
      herm1(1,1) = cplx(2,0);
      herm1(2,2) = cplx(2,0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 *= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 2, 0) || herm2(0,1) != cplx(-8,-2) || herm2(0,2) != cplx(14,6) ||
          herm2(1,0) != cplx(-8, 2) || herm2(1,1) != cplx( 4, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(14,-6) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix multiplication assignment (HermitianMatrix)";

      OHT herm1( 3UL, 3UL );
      herm1(0,0) = cplx(2,0);
      herm1(1,1) = cplx(2,0);
      herm1(2,2) = cplx(2,0);

      HT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 *= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 2, 0) || herm2(0,1) != cplx(-8,-2) || herm2(0,2) != cplx(14,6) ||
          herm2(1,0) != cplx(-8, 2) || herm2(1,1) != cplx( 4, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(14,-6) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major dense matrix multiplication assignment (Hermitian)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix multiplication assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,0) = cplx(2,0);
      mat(1,1) = cplx(2,0);
      mat(2,2) = cplx(2,0);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm *= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 2, 0) || herm(0,1) != cplx(-8,-2) || herm(0,2) != cplx(14,6) ||
          herm(1,0) != cplx(-8, 2) || herm(1,1) != cplx( 4, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(14,-6) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (Hermitian)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix multiplication assignment (Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,0) = cplx(2,0);
      mat(1,1) = cplx(2,0);
      mat(2,2) = cplx(2,0);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm *= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 2, 0) || herm(0,1) != cplx(-8,-2) || herm(0,2) != cplx(14,6) ||
          herm(1,0) != cplx(-8, 2) || herm(1,1) != cplx( 4, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(14,-6) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix multiplication assignment (non-Hermitian)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix multiplication assignment (non-Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-2);
      mat(0,2) = cplx( 6, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 3);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1,0);
      herm(0,1) = cplx(-4,0);
      herm(0,2) = cplx( 7,0);
      herm(1,1) = cplx( 2,0);
      herm(2,2) = cplx( 3,0);

      try {
         herm *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix multiplication assignment (non-Hermitian)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix multiplication assignment (non-Hermitian)";

      blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, cplx(0,0) );
      mat(0,1) = cplx(-2,-2);
      mat(0,2) = cplx( 6, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 3);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1,0);
      herm(0,1) = cplx(-4,0);
      herm(0,2) = cplx( 7,0);
      herm(1,1) = cplx( 2,0);
      herm(2,2) = cplx( 3,0);

      try {
         herm *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix multiplication assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix multiplication assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::rowMajor> > herm1( 3UL );
      herm1(0,0) = cplx(2,0);
      herm1(1,1) = cplx(2,0);
      herm1(2,2) = cplx(2,0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 *= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 2, 0) || herm2(0,1) != cplx(-8,-2) || herm2(0,2) != cplx(14,6) ||
          herm2(1,0) != cplx(-8, 2) || herm2(1,1) != cplx( 4, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(14,-6) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix multiplication assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::columnMajor> > herm1( 3UL );
      herm1(0,0) = cplx(2,0);
      herm1(1,1) = cplx(2,0);
      herm1(2,2) = cplx(2,0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 *= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 2, 0) || herm2(0,1) != cplx(-8,-2) || herm2(0,2) != cplx(14,6) ||
          herm2(1,0) != cplx(-8, 2) || herm2(1,1) != cplx( 4, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(14,-6) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix multiplication assignment (Hermitian)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix multiplication assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = cplx(2,0);
      mat(1,1) = cplx(2,0);
      mat(2,2) = cplx(2,0);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm *= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 2, 0) || herm(0,1) != cplx(-8,-2) || herm(0,2) != cplx(14,6) ||
          herm(1,0) != cplx(-8, 2) || herm(1,1) != cplx( 4, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(14,-6) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (Hermitian)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix multiplication assignment (Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = cplx(2,0);
      mat(1,1) = cplx(2,0);
      mat(2,2) = cplx(2,0);
      mat.insert( 1UL, 2UL, cplx(0,0) );

      OHT herm( 3UL );
      herm(0,0) = cplx( 1, 0);
      herm(0,1) = cplx(-4,-1);
      herm(0,2) = cplx( 7, 3);
      herm(1,1) = cplx( 2, 0);
      herm(2,2) = cplx( 3, 0);

      herm *= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 7UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != cplx( 2, 0) || herm(0,1) != cplx(-8,-2) || herm(0,2) != cplx(14,6) ||
          herm(1,0) != cplx(-8, 2) || herm(1,1) != cplx( 4, 0) || herm(1,2) != cplx( 0,0) ||
          herm(2,0) != cplx(14,-6) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix multiplication assignment (non-Hermitian)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix multiplication assignment (non-Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = cplx(-2,-2);
      mat(0,2) = cplx( 6, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 3);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1,0);
      herm(0,1) = cplx(-4,0);
      herm(0,2) = cplx( 7,0);
      herm(1,1) = cplx( 2,0);
      herm(2,2) = cplx( 3,0);

      try {
         herm *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix multiplication assignment (non-Hermitian)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix multiplication assignment (non-Hermitian)";

      blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = cplx(-2,-2);
      mat(0,2) = cplx( 6, 1);
      mat(1,1) = cplx( 3, 0);
      mat(2,0) = cplx( 6, 3);

      OHT herm( 3UL );
      herm(0,0) = cplx( 1,0);
      herm(0,1) = cplx(-4,0);
      herm(0,2) = cplx( 7,0);
      herm(1,1) = cplx( 2,0);
      herm(2,2) = cplx( 3,0);

      try {
         herm *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-Hermitian row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix multiplication assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix multiplication assignment (HermitianMatrix)";

      HT herm1( 3UL, 3UL );
      herm1(0,0) = cplx(2,0);
      herm1(1,1) = cplx(2,0);
      herm1(2,2) = cplx(2,0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 *= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 2, 0) || herm2(0,1) != cplx(-8,-2) || herm2(0,2) != cplx(14,6) ||
          herm2(1,0) != cplx(-8, 2) || herm2(1,1) != cplx( 4, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(14,-6) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix multiplication assignment (HermitianMatrix)";

      OHT herm1( 3UL, 3UL );
      herm1(0,0) = cplx(2,0);
      herm1(1,1) = cplx(2,0);
      herm1(2,2) = cplx(2,0);

      OHT herm2( 3UL );
      herm2(0,0) = cplx( 1, 0);
      herm2(0,1) = cplx(-4,-1);
      herm2(0,2) = cplx( 7, 3);
      herm2(1,1) = cplx( 2, 0);
      herm2(2,2) = cplx( 3, 0);

      herm2 *= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 7UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) != cplx( 2, 0) || herm2(0,1) != cplx(-8,-2) || herm2(0,2) != cplx(14,6) ||
          herm2(1,0) != cplx(-8, 2) || herm2(1,1) != cplx( 4, 0) || herm2(1,2) != cplx( 0,0) ||
          herm2(2,0) != cplx(14,-6) || herm2(2,1) != cplx( 0, 0) || herm2(2,2) != cplx( 6,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n"
                                     "( -8  4  0 )\n"
                                     "( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace hermitianmatrix

} // namespace adaptors

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
   std::cout << "   Running HermitianMatrix sparse complex test (part 1)..." << std::endl;

   try
   {
      RUN_HERMITIANMATRIX_SPARSECOMPLEX_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during HermitianMatrix sparse complex test (part 1):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
