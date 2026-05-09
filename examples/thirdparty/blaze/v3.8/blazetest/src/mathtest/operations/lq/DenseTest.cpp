//=================================================================================================
/*!
//  \file src/mathtest/operations/lq/DenseTest.cpp
//  \brief Source file for the dense matrix LQ test
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
#include <blazetest/mathtest/operations/lq/DenseTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace operations {

namespace lq {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DenseTest test.
//
// \exception std::runtime_error Error during LQ decomposition detected.
*/
DenseTest::DenseTest()
{
   testGeneral();
   testSymmetric();
   testHermitian();
   testLower();
   testUniLower();
   testUpper();
   testUniUpper();
   testDiagonal();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the LQ decomposition functionality for general matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LQ decomposition for general matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testGeneral()
{
   //testRandom< blaze::DynamicMatrix<float,blaze::rowMajor   > >();
   //testRandom< blaze::DynamicMatrix<float,blaze::columnMajor> >();

   testRandom< blaze::DynamicMatrix<double,blaze::rowMajor   > >();
   testRandom< blaze::DynamicMatrix<double,blaze::columnMajor> >();

   //testRandom< blaze::DynamicMatrix<cfloat,blaze::rowMajor   > >();
   //testRandom< blaze::DynamicMatrix<cfloat,blaze::columnMajor> >();

   testRandom< blaze::DynamicMatrix<cdouble,blaze::rowMajor   > >();
   testRandom< blaze::DynamicMatrix<cdouble,blaze::columnMajor> >();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LQ decomposition functionality for symmetric matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LQ decomposition for symmetric matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testSymmetric()
{
   //testRandom< blaze::SymmetricMatrix< blaze::DynamicMatrix<float,blaze::rowMajor   > > >();
   //testRandom< blaze::SymmetricMatrix< blaze::DynamicMatrix<float,blaze::columnMajor> > >();

   testRandom< blaze::SymmetricMatrix< blaze::DynamicMatrix<double,blaze::rowMajor   > > >();
   testRandom< blaze::SymmetricMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > >();

   //testRandom< blaze::SymmetricMatrix< blaze::DynamicMatrix<cfloat,blaze::rowMajor   > > >();
   //testRandom< blaze::SymmetricMatrix< blaze::DynamicMatrix<cfloat,blaze::columnMajor> > >();

   testRandom< blaze::SymmetricMatrix< blaze::DynamicMatrix<cdouble,blaze::rowMajor   > > >();
   testRandom< blaze::SymmetricMatrix< blaze::DynamicMatrix<cdouble,blaze::columnMajor> > >();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LQ decomposition functionality for Hermitian matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LQ decomposition for Hermitian matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testHermitian()
{
   //testRandom< blaze::HermitianMatrix< blaze::DynamicMatrix<float,blaze::rowMajor   > > >();
   //testRandom< blaze::HermitianMatrix< blaze::DynamicMatrix<float,blaze::columnMajor> > >();

   testRandom< blaze::HermitianMatrix< blaze::DynamicMatrix<double,blaze::rowMajor   > > >();
   testRandom< blaze::HermitianMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > >();

   //testRandom< blaze::HermitianMatrix< blaze::DynamicMatrix<cfloat,blaze::rowMajor   > > >();
   //testRandom< blaze::HermitianMatrix< blaze::DynamicMatrix<cfloat,blaze::columnMajor> > >();

   testRandom< blaze::HermitianMatrix< blaze::DynamicMatrix<cdouble,blaze::rowMajor   > > >();
   testRandom< blaze::HermitianMatrix< blaze::DynamicMatrix<cdouble,blaze::columnMajor> > >();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LQ decomposition functionality for lower matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LQ decomposition for lower matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testLower()
{
   //testRandom< blaze::LowerMatrix< blaze::DynamicMatrix<float,blaze::rowMajor   > > >();
   //testRandom< blaze::LowerMatrix< blaze::DynamicMatrix<float,blaze::columnMajor> > >();

   testRandom< blaze::LowerMatrix< blaze::DynamicMatrix<double,blaze::rowMajor   > > >();
   testRandom< blaze::LowerMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > >();

   //testRandom< blaze::LowerMatrix< blaze::DynamicMatrix<cfloat,blaze::rowMajor   > > >();
   //testRandom< blaze::LowerMatrix< blaze::DynamicMatrix<cfloat,blaze::columnMajor> > >();

   testRandom< blaze::LowerMatrix< blaze::DynamicMatrix<cdouble,blaze::rowMajor   > > >();
   testRandom< blaze::LowerMatrix< blaze::DynamicMatrix<cdouble,blaze::columnMajor> > >();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LQ decomposition functionality for unilower matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LQ decomposition for unilower matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testUniLower()
{
   //testRandom< blaze::UniLowerMatrix< blaze::DynamicMatrix<float,blaze::rowMajor   > > >();
   //testRandom< blaze::UniLowerMatrix< blaze::DynamicMatrix<float,blaze::columnMajor> > >();

   testRandom< blaze::UniLowerMatrix< blaze::DynamicMatrix<double,blaze::rowMajor   > > >();
   testRandom< blaze::UniLowerMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > >();

   //testRandom< blaze::UniLowerMatrix< blaze::DynamicMatrix<cfloat,blaze::rowMajor   > > >();
   //testRandom< blaze::UniLowerMatrix< blaze::DynamicMatrix<cfloat,blaze::columnMajor> > >();

   testRandom< blaze::UniLowerMatrix< blaze::DynamicMatrix<cdouble,blaze::rowMajor   > > >();
   testRandom< blaze::UniLowerMatrix< blaze::DynamicMatrix<cdouble,blaze::columnMajor> > >();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LQ decomposition functionality for upper matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LQ decomposition for upper matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testUpper()
{
   //testRandom< blaze::UpperMatrix< blaze::DynamicMatrix<float,blaze::rowMajor   > > >();
   //testRandom< blaze::UpperMatrix< blaze::DynamicMatrix<float,blaze::columnMajor> > >();

   testRandom< blaze::UpperMatrix< blaze::DynamicMatrix<double,blaze::rowMajor   > > >();
   testRandom< blaze::UpperMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > >();

   //testRandom< blaze::UpperMatrix< blaze::DynamicMatrix<cfloat,blaze::rowMajor   > > >();
   //testRandom< blaze::UpperMatrix< blaze::DynamicMatrix<cfloat,blaze::columnMajor> > >();

   testRandom< blaze::UpperMatrix< blaze::DynamicMatrix<cdouble,blaze::rowMajor   > > >();
   testRandom< blaze::UpperMatrix< blaze::DynamicMatrix<cdouble,blaze::columnMajor> > >();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LQ decomposition functionality for uniupper matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LQ decomposition for uniupper matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testUniUpper()
{
   //testRandom< blaze::UniUpperMatrix< blaze::DynamicMatrix<float,blaze::rowMajor   > > >();
   //testRandom< blaze::UniUpperMatrix< blaze::DynamicMatrix<float,blaze::columnMajor> > >();

   testRandom< blaze::UniUpperMatrix< blaze::DynamicMatrix<double,blaze::rowMajor   > > >();
   testRandom< blaze::UniUpperMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > >();

   //testRandom< blaze::UniUpperMatrix< blaze::DynamicMatrix<cfloat,blaze::rowMajor   > > >();
   //testRandom< blaze::UniUpperMatrix< blaze::DynamicMatrix<cfloat,blaze::columnMajor> > >();

   testRandom< blaze::UniUpperMatrix< blaze::DynamicMatrix<cdouble,blaze::rowMajor   > > >();
   testRandom< blaze::UniUpperMatrix< blaze::DynamicMatrix<cdouble,blaze::columnMajor> > >();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the LQ decomposition functionality for diagonal matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LQ decomposition for diagonal matrices. In case an error
// is detected, a \a std::runtime_error exception is thrown.
*/
void DenseTest::testDiagonal()
{
   //testRandom< blaze::DiagonalMatrix< blaze::DynamicMatrix<float,blaze::rowMajor   > > >();
   //testRandom< blaze::DiagonalMatrix< blaze::DynamicMatrix<float,blaze::columnMajor> > >();

   testRandom< blaze::DiagonalMatrix< blaze::DynamicMatrix<double,blaze::rowMajor   > > >();
   testRandom< blaze::DiagonalMatrix< blaze::DynamicMatrix<double,blaze::columnMajor> > >();

   //testRandom< blaze::DiagonalMatrix< blaze::DynamicMatrix<cfloat,blaze::rowMajor   > > >();
   //testRandom< blaze::DiagonalMatrix< blaze::DynamicMatrix<cfloat,blaze::columnMajor> > >();

   testRandom< blaze::DiagonalMatrix< blaze::DynamicMatrix<cdouble,blaze::rowMajor   > > >();
   testRandom< blaze::DiagonalMatrix< blaze::DynamicMatrix<cdouble,blaze::columnMajor> > >();
}
//*************************************************************************************************

} // namespace lq

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
   std::cout << "   Running dense matrix LQ decomposition test..." << std::endl;

   try
   {
      RUN_DENSE_LQ_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during dense matrix LQ decomposition test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
