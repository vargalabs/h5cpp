//=================================================================================================
/*!
//  \file src/mathtest/operations/svecsvecmult/AliasingTest.cpp
//  \brief Source file for the sparse vector/sparse vector multiplication aliasing test
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
#include <blazetest/mathtest/operations/svecsvecmult/AliasingTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace operations {

namespace svecsvecmult {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the aliasing test class.
//
// \exception std::runtime_error Operation error detected.
*/
AliasingTest::AliasingTest()
   : sa4_   ( 4UL )
   , sb3_   ( 3UL )
   , sc3_   ( 3UL )
   , sA3x4_ ( 3UL, 4UL )
   , sB3x3_ ( 3UL, 3UL )
   , result_()
   , test_  ()
{
   testSVecSVecMult();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the sparse vector/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs aliasing tests for the sparse vector/sparse vector multiplication.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AliasingTest::testSVecSVecMult()
{
   //=====================================================================================
   // Multiplication
   //=====================================================================================

   // Assignment to left-hand side operand (1)
   {
      test_ = "SVecSVecMult - Assignment to left-hand side operand (1)";

      initialize();

      result_ = sb3_ * sc3_;
      sb3_    = sb3_ * sc3_;

      checkResult( sb3_, result_ );
   }

   // Assignment to left-hand side operand (2)
   {
      test_ = "SVecSVecMult - Assignment to left-hand side operand (2)";

      initialize();

      result_ = sb3_ * eval( sc3_ );
      sb3_    = sb3_ * eval( sc3_ );

      checkResult( sb3_, result_ );
   }

   // Assignment to left-hand side compound operand
   {
      test_ = "SVecSVecMult - Assignment to left-hand side compound operand";

      initialize();

      result_ = ( sA3x4_ * sa4_ ) * sc3_;
      sa4_    = ( sA3x4_ * sa4_ ) * sc3_;

      checkResult( sa4_, result_ );
   }

   // Assignment to right-hand side operand (1)
   {
      test_ = "SVecSVecMult - Assignment to right-hand side operand (1)";

      initialize();

      result_ = sb3_ * sc3_;
      sc3_    = sb3_ * sc3_;

      checkResult( sc3_, result_ );
   }

   // Assignment to right-hand side operand (2)
   {
      test_ = "SVecSVecMult - Assignment to right-hand side operand (2)";

      initialize();

      result_ = eval( sb3_ ) * sc3_;
      sc3_    = eval( sb3_ ) * sc3_;

      checkResult( sc3_, result_ );
   }

   // Assignment to right-hand side compound operand
   {
      test_ = "SVecSVecMult - Assignment to right-hand side compound operand";

      initialize();

      result_ = sb3_ * ( sA3x4_ * sa4_ );
      sa4_    = sb3_ * ( sA3x4_ * sa4_ );

      checkResult( sa4_, result_ );
   }

   // Complex operation: a = ( 2*a ) * ( A * b );
   {
      test_ = "SVecSVecMult - Complex operation: a = ( 2*a ) * ( A * b );";

      initialize();

      result_ = ( 2*sb3_ ) * ( sA3x4_ * sa4_ );
      sb3_    = ( 2*sb3_ ) * ( sA3x4_ * sa4_ );

      checkResult( sb3_, result_ );
   }

   // Complex operation: a = ( A * b ) * ( 2*a );
   {
      test_ = "SVecSVecMult - Complex operation: a = ( A * b ) * ( 2*a );";

      initialize();

      result_ = ( sA3x4_ * sa4_ ) * ( 2*sb3_ );
      sb3_    = ( sA3x4_ * sa4_ ) * ( 2*sb3_ );

      checkResult( sb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with addition assignment
   //=====================================================================================

   // Addition assignment to left-hand side operand (1)
   {
      test_ = "SVecSVecMult - Addition assignment to left-hand side operand (1)";

      initialize();

      result_ =  sb3_;
      result_ += sb3_ * sc3_;
      sb3_    += sb3_ * sc3_;

      checkResult( sb3_, result_ );
   }

   // Addition assignment to left-hand side operand (2)
   {
      test_ = "SVecSVecMult - Addition assignment to left-hand side operand (2)";

      initialize();

      result_ =  sb3_;
      result_ += sb3_ * eval( sc3_ );
      sb3_    += sb3_ * eval( sc3_ );

      checkResult( sb3_, result_ );
   }

   // Addition assignment to left-hand side dense compound operand
   {
      test_ = "SVecSVecMult - Addition assignment to left-hand side dense compound operand";

      initialize();

      result_ =  sb3_;
      result_ += ( sB3x3_ * sb3_ ) * sc3_;
      sb3_    += ( sB3x3_ * sb3_ ) * sc3_;

      checkResult( sb3_, result_ );
   }

   // Addition assignment to right-hand side operand (1)
   {
      test_ = "SVecSVecMult - Addition assignment to right-hand side operand (1)";

      initialize();

      result_ =  sc3_;
      result_ += sb3_ * sc3_;
      sc3_    += sb3_ * sc3_;

      checkResult( sc3_, result_ );
   }

   // Addition assignment to right-hand side operand (2)
   {
      test_ = "SVecSVecMult - Addition assignment to right-hand side operand (2)";

      initialize();

      result_ =  sc3_;
      result_ += eval( sb3_ ) * sc3_;
      sc3_    += eval( sb3_ ) * sc3_;

      checkResult( sc3_, result_ );
   }

   // Addition assignment to right-hand side dense compound operand
   {
      test_ = "SVecSVecMult - Addition assignment to right-hand side dense compound operand";

      initialize();

      result_ =  sc3_;
      result_ += sb3_ * ( sB3x3_ * sc3_ );
      sc3_    += sb3_ * ( sB3x3_ * sc3_ );

      checkResult( sc3_, result_ );
   }

   // Complex operation: a += ( 2*a ) * ( A * b );
   {
      test_ = "SVecSVecMult - Complex operation: a += ( 2*a ) * ( A * b );";

      initialize();

      result_ =  sb3_;
      result_ += ( 2*sb3_ ) * ( sA3x4_ * sa4_ );
      sb3_    += ( 2*sb3_ ) * ( sA3x4_ * sa4_ );

      checkResult( sb3_, result_ );
   }

   // Complex operation: a += ( A * b ) * ( 2*a );
   {
      test_ = "SVecSVecMult - Complex operation: a += ( A * b ) * ( 2*a );";

      initialize();

      result_ =  sb3_;
      result_ += ( sA3x4_ * sa4_ ) * ( 2*sb3_ );
      sb3_    += ( sA3x4_ * sa4_ ) * ( 2*sb3_ );

      checkResult( sb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with subtraction assignment
   //=====================================================================================

   // Subtraction assignment to left-hand side operand (1)
   {
      test_ = "SVecSVecMult - Subtraction assignment to left-hand side operand (1)";

      initialize();

      result_ =  sb3_;
      result_ -= sb3_ * sc3_;
      sb3_    -= sb3_ * sc3_;

      checkResult( sb3_, result_ );
   }

   // Subtraction assignment to left-hand side operand (2)
   {
      test_ = "SVecSVecMult - Subtraction assignment to left-hand side operand (2)";

      initialize();

      result_ =  sb3_;
      result_ -= sb3_ * eval( sc3_ );
      sb3_    -= sb3_ * eval( sc3_ );

      checkResult( sb3_, result_ );
   }

   // Subtraction assignment to left-hand side dense compound operand
   {
      test_ = "SVecSVecMult - Subtraction assignment to left-hand side dense compound operand";

      initialize();

      result_ =  sb3_;
      result_ -= ( sB3x3_ * sb3_ ) * sc3_;
      sb3_    -= ( sB3x3_ * sb3_ ) * sc3_;

      checkResult( sb3_, result_ );
   }

   // Subtraction assignment to right-hand side operand (1)
   {
      test_ = "SVecSVecMult - Subtraction assignment to right-hand side operand (1)";

      initialize();

      result_ =  sc3_;
      result_ -= sb3_ * sc3_;
      sc3_    -= sb3_ * sc3_;

      checkResult( sc3_, result_ );
   }

   // Subtraction assignment to right-hand side operand (2)
   {
      test_ = "SVecSVecMult - Subtraction assignment to right-hand side operand (2)";

      initialize();

      result_ =  sc3_;
      result_ -= eval( sb3_ ) * sc3_;
      sc3_    -= eval( sb3_ ) * sc3_;

      checkResult( sc3_, result_ );
   }

   // Subtraction assignment to right-hand side dense compound operand
   {
      test_ = "SVecSVecMult - Subtraction assignment to right-hand side dense compound operand";

      initialize();

      result_ =  sc3_;
      result_ -= sb3_ * ( sB3x3_ * sc3_ );
      sc3_    -= sb3_ * ( sB3x3_ * sc3_ );

      checkResult( sc3_, result_ );
   }

   // Complex operation: a -= ( 2*a ) * ( A * b );
   {
      test_ = "SVecSVecMult - Complex operation: a -= ( 2*a ) * ( A * b );";

      initialize();

      result_ =  sb3_;
      result_ -= ( 2*sb3_ ) * ( sA3x4_ * sa4_ );
      sb3_    -= ( 2*sb3_ ) * ( sA3x4_ * sa4_ );

      checkResult( sb3_, result_ );
   }

   // Complex operation: a -= ( A * b ) * ( 2*a );
   {
      test_ = "SVecSVecMult - Complex operation: a -= ( A * b ) * ( 2*a );";

      initialize();

      result_ =  sb3_;
      result_ -= ( sA3x4_ * sa4_ ) * ( 2*sb3_ );
      sb3_    -= ( sA3x4_ * sa4_ ) * ( 2*sb3_ );

      checkResult( sb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with multiplication assignment
   //=====================================================================================

   // Multiplication assignment to left-hand side operand (1)
   {
      test_ = "SVecSVecMult - Multiplication assignment to left-hand side operand (1)";

      initialize();

      result_ =  sb3_;
      result_ *= sb3_ * sc3_;
      sb3_    *= sb3_ * sc3_;

      checkResult( sb3_, result_ );
   }

   // Multiplication assignment to left-hand side operand (2)
   {
      test_ = "SVecSVecMult - Multiplication assignment to left-hand side operand (2)";

      initialize();

      result_ =  sb3_;
      result_ *= sb3_ * eval( sc3_ );
      sb3_    *= sb3_ * eval( sc3_ );

      checkResult( sb3_, result_ );
   }

   // Multiplication assignment to left-hand side dense compound operand
   {
      test_ = "SVecSVecMult - Multiplication assignment to left-hand side dense compound operand";

      initialize();

      result_ =  sb3_;
      result_ *= ( sB3x3_ * sb3_ ) * sc3_;
      sb3_    *= ( sB3x3_ * sb3_ ) * sc3_;

      checkResult( sb3_, result_ );
   }

   // Multiplication assignment to right-hand side operand (1)
   {
      test_ = "SVecSVecMult - Multiplication assignment to right-hand side operand (1)";

      initialize();

      result_ =  sc3_;
      result_ *= sb3_ * sc3_;
      sc3_    *= sb3_ * sc3_;

      checkResult( sc3_, result_ );
   }

   // Multiplication assignment to right-hand side operand (2)
   {
      test_ = "SVecSVecMult - Multiplication assignment to right-hand side operand (2)";

      initialize();

      result_ =  sc3_;
      result_ *= eval( sb3_ ) * sc3_;
      sc3_    *= eval( sb3_ ) * sc3_;

      checkResult( sc3_, result_ );
   }

   // Multiplication assignment to right-hand side dense compound operand
   {
      test_ = "SVecSVecMult - Multiplication assignment to right-hand side dense compound operand";

      initialize();

      result_ =  sc3_;
      result_ *= sb3_ * ( sB3x3_ * sc3_ );
      sc3_    *= sb3_ * ( sB3x3_ * sc3_ );

      checkResult( sc3_, result_ );
   }

   // Complex operation: a *= ( 2*a ) * ( A * b );
   {
      test_ = "SVecSVecMult - Complex operation: a *= ( 2*a ) * ( A * b );";

      initialize();

      result_ =  sb3_;
      result_ *= ( 2*sb3_ ) * ( sA3x4_ * sa4_ );
      sb3_    *= ( 2*sb3_ ) * ( sA3x4_ * sa4_ );

      checkResult( sb3_, result_ );
   }

   // Complex operation: a *= ( A * b ) * ( 2*a );
   {
      test_ = "SVecSVecMult - Complex operation: a *= ( A * b ) * ( 2*a );";

      initialize();

      result_ =  sb3_;
      result_ *= ( sA3x4_ * sa4_ ) * ( 2*sb3_ );
      sb3_    *= ( sA3x4_ * sa4_ ) * ( 2*sb3_ );

      checkResult( sb3_, result_ );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization of all member vectors and matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function initializes all member vectors and matrices to specific predetermined values.
*/
void AliasingTest::initialize()
{
   //=====================================================================================
   // Initialization of the sparse vectors
   //=====================================================================================

   // Initializing the first sparse column vector
   sa4_.resize( 4UL, false );
   sa4_.reset();
   sa4_[0] = -1;
   sa4_[2] = -3;
   sa4_[3] =  2;

   // Initializing the second sparse column vector
   sb3_.resize( 3UL, false );
   sb3_.reset();
   sb3_[0] = 1;
   sb3_[1] = 2;
   sb3_[2] = 3;

   // Initializing the third sparse column vector
   sc3_.resize( 3UL, false );
   sc3_.reset();
   sc3_[1] = 2;
   sc3_[2] = 1;


   //=====================================================================================
   // Initialization of the sparse matrices
   //=====================================================================================

   // Initializing the first row-major sparse matrix
   sA3x4_(0,0) = -1;
   sA3x4_(0,2) = -2;
   sA3x4_(1,1) =  2;
   sA3x4_(1,2) = -3;
   sA3x4_(1,3) =  1;
   sA3x4_(2,1) =  1;
   sA3x4_(2,2) =  2;
   sA3x4_(2,3) =  2;

   // Initializing the second row-major sparse matrix
   sB3x3_(0,0) = -1;
   sB3x3_(1,0) =  1;
   sB3x3_(1,1) = -2;
   sB3x3_(1,2) =  2;
   sB3x3_(2,2) = -3;
}
//*************************************************************************************************

} // namespace svecsvecmult

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
   std::cout << "   Running aliasing test..." << std::endl;

   try
   {
      RUN_SVECSVECMULT_ALIASING_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during aliasing test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
