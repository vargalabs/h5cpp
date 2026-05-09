//=================================================================================================
/*!
//  \file src/mathtest/operations/inversion/DenseTest.cpp
//  \brief Source file for the dense matrix inversion test
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
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/operations/inversion/DenseTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace operations {

namespace inversion {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DenseTest inversion test.
//
// \exception std::runtime_error Matrix inversion error detected.
*/
DenseTest::DenseTest()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::HermitianMatrix;
   using blaze::LowerMatrix;
   using blaze::UniLowerMatrix;
   using blaze::UpperMatrix;
   using blaze::UniUpperMatrix;
   using blaze::DiagonalMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   using cdouble = blaze::complex<double>;


   //=====================================================================================
   // Specific matrix tests
   //=====================================================================================

   testSpecific();


   //=====================================================================================
   // Random matrix tests
   //=====================================================================================

   for( size_t i=0UL; i<12UL; ++i )
   {
      testRandom< DynamicMatrix<double ,rowMajor> >( i );
      testRandom< DynamicMatrix<cdouble,rowMajor> >( i );

      testRandom< DynamicMatrix<double ,columnMajor> >( i );
      testRandom< DynamicMatrix<cdouble,columnMajor> >( i );

      testRandom< SymmetricMatrix< DynamicMatrix<double ,rowMajor> > >( i );
      testRandom< SymmetricMatrix< DynamicMatrix<cdouble,rowMajor> > >( i );
      testRandom< HermitianMatrix< DynamicMatrix<double ,rowMajor> > >( i );
      testRandom< HermitianMatrix< DynamicMatrix<cdouble,rowMajor> > >( i );
      testRandom< LowerMatrix    < DynamicMatrix<double ,rowMajor> > >( i );
      testRandom< UniLowerMatrix < DynamicMatrix<double ,rowMajor> > >( i );
      testRandom< UpperMatrix    < DynamicMatrix<double ,rowMajor> > >( i );
      testRandom< UniUpperMatrix < DynamicMatrix<double ,rowMajor> > >( i );
      testRandom< DiagonalMatrix < DynamicMatrix<double ,rowMajor> > >( i );

      testRandom< SymmetricMatrix< DynamicMatrix<double ,columnMajor> > >( i );
      testRandom< SymmetricMatrix< DynamicMatrix<cdouble,columnMajor> > >( i );
      testRandom< HermitianMatrix< DynamicMatrix<double ,columnMajor> > >( i );
      testRandom< HermitianMatrix< DynamicMatrix<cdouble,columnMajor> > >( i );
      testRandom< LowerMatrix    < DynamicMatrix<double ,columnMajor> > >( i );
      testRandom< UniLowerMatrix < DynamicMatrix<double ,columnMajor> > >( i );
      testRandom< UpperMatrix    < DynamicMatrix<double ,columnMajor> > >( i );
      testRandom< UniUpperMatrix < DynamicMatrix<double ,columnMajor> > >( i );
      testRandom< DiagonalMatrix < DynamicMatrix<double ,columnMajor> > >( i );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the determinant functionality with specific, predetermined matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function computes determinants for specific, predetermined matrices. In case an error is
// detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSpecific()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   //=====================================================================================
   // Row-major matrix tests
   //=====================================================================================

   {
      {
         test_ = "Row-major dense matrix inversion (0x0)";

         blaze::DynamicMatrix<double,blaze::rowMajor> A;

         invert( A );

         if( A.rows() != 0UL || A.columns() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix inversion failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Row-major dense matrix inversion (non-square)";

         blaze::DynamicMatrix<double,blaze::rowMajor> A( 2UL, 3UL );

         try {
            invert( A );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inversion of a non-square matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major matrix tests
   //=====================================================================================

   {
      {
         test_ = "Column-major dense matrix inversion (0x0)";

         blaze::DynamicMatrix<double,blaze::columnMajor> A;

         invert( A );

         if( A.rows() != 0UL || A.columns() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Matrix inversion failed\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         test_ = "Column-major dense matrix inversion (non-square)";

         blaze::DynamicMatrix<double,blaze::columnMajor> A( 2UL, 3UL );

         try {
            invert( A );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inversion of a non-square matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << A << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

#endif
}
//*************************************************************************************************

} // namespace inversion

} // namespace operations

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
   std::cout << "   Running dense matrix inversion test..." << std::endl;

   try
   {
      RUN_INVERSION_DENSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix inversion test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
