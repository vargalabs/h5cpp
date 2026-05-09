//=================================================================================================
/*!
//  \file src/mathtest/operations/svecdvecsub/AliasingTest.cpp
//  \brief Source file for the dense vector/sparse vector subtraction aliasing test
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
#include <blaze/math/Subvector.h>
#include <blazetest/mathtest/operations/svecdvecsub/AliasingTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace operations {

namespace svecdvecsub {

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
   , da4_   ( 4UL )
   , db3_   ( 3UL )
   , sA3x4_ ( 3UL, 4UL )
   , sB3x3_ ( 3UL, 3UL )
   , dA3x4_ ( 3UL, 4UL )
   , dB3x3_ ( 3UL, 3UL )
   , result_()
   , test_  ()
{
   testSVecDVecSub();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the dense vector/sparse vector subtraction.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs aliasing tests for the dense vector/sparse vector subtraction.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AliasingTest::testSVecDVecSub()
{
   using blaze::subvector;


   //=====================================================================================
   // Subtraction
   //=====================================================================================

   // Assignment to left-hand side operand (1)
   {
      test_ = "SVecDVecSub - Assignment to left-hand side operand (1)";

      initialize();

      result_ = sb3_ - db3_;
      sb3_    = sb3_ - db3_;

      checkResult( sb3_, result_ );
   }

   // Assignment to left-hand side operand (2)
   {
      test_ = "SVecDVecSub - Assignment to left-hand side operand (2)";

      initialize();

      result_ = sb3_ - eval( db3_ );
      sb3_    = sb3_ - eval( db3_ );

      checkResult( sb3_, result_ );
   }

   // Assignment to left-hand side compound operand
   {
      test_ = "SVecDVecSub - Assignment to left-hand side compound operand";

      initialize();

      result_ = ( sA3x4_ * sa4_ ) - db3_;
      sa4_   =  ( sA3x4_ * sa4_ ) - db3_;

      checkResult( sa4_, result_ );
   }

   // Assignment to left-hand side subvector operand
   {
      test_ = "SVecDVecSub - Assignment to left-hand side subvector operand";

      initialize();

      result_ = subvector( sa4_, 1UL, 3UL ) - db3_;
      sa4_    = subvector( sa4_, 1UL, 3UL ) - db3_;

      checkResult( sa4_, result_ );
   }

   // Assignment to right-hand side operand (1)
   {
      test_ = "SVecDVecSub - Assignment to right-hand side operand (1)";

      initialize();

      result_ = sb3_ - db3_;
      db3_    = sb3_ - db3_;

      checkResult( db3_, result_ );
   }

   // Assignment to right-hand side operand (2)
   {
      test_ = "SVecDVecSub - Assignment to right-hand side operand (2)";

      initialize();

      result_ = eval( sb3_ ) - db3_;
      db3_    = eval( sb3_ ) - db3_;

      checkResult( db3_, result_ );
   }

   // Assignment to right-hand side compound operand
   {
      test_ = "SVecDVecSub - Assignment to right-hand side compound operand";

      initialize();

      result_ = sb3_ - ( dA3x4_ * da4_ );
      da4_    = sb3_ - ( dA3x4_ * da4_ );

      checkResult( da4_, result_ );
   }

   // Assignment to right-hand side subvector operand
   {
      test_ = "SVecDVecSub - Assignment to right-hand side subvector operand";

      initialize();

      result_ = sb3_ - subvector( da4_, 1UL, 3UL );
      da4_    = sb3_ - subvector( da4_, 1UL, 3UL );

      checkResult( da4_, result_ );
   }

   // Complex operation: a = ( 2*a ) - ( A * b );
   {
      test_ = "SVecDVecSub - Complex operation: a = ( 2*a ) - ( A * b );";

      initialize();

      result_ = ( 2*sb3_ ) - ( dA3x4_ * da4_ );
      sb3_    = ( 2*sb3_ ) - ( dA3x4_ * da4_ );

      checkResult( sb3_, result_ );
   }

   // Complex operation: a = ( A * b ) - ( 2*a );
   {
      test_ = "SVecDVecSub - Complex operation: a = ( A * b ) - ( 2*a );";

      initialize();

      result_ = ( sA3x4_ * sa4_ ) - ( 2*db3_ );
      db3_    = ( sA3x4_ * sa4_ ) - ( 2*db3_ );

      checkResult( db3_, result_ );
   }


   //=====================================================================================
   // Subtraction with addition assignment
   //=====================================================================================

   // Addition assignment to left-hand side operand (1)
   {
      test_ = "SVecDVecSub - Addition assignment to left-hand side operand (1)";

      initialize();

      result_ =  sb3_;
      result_ += sb3_ - db3_;
      sb3_    += sb3_ - db3_;

      checkResult( sb3_, result_ );
   }

   // Addition assignment to left-hand side operand (2)
   {
      test_ = "SVecDVecSub - Addition assignment to left-hand side operand (2)";

      initialize();

      result_ =  sb3_;
      result_ += sb3_ - eval( db3_ );
      sb3_    += sb3_ - eval( db3_ );

      checkResult( sb3_, result_ );
   }

   // Addition assignment to left-hand side compound operand
   {
      test_ = "SVecDVecSub - Addition assignment to left-hand side compound operand";

      initialize();

      result_ =  sb3_;
      result_ += ( sB3x3_ * sb3_ ) - db3_;
      sb3_    += ( sB3x3_ * sb3_ ) - db3_;

      checkResult( sb3_, result_ );
   }

   // Addition assignment to right-hand side operand (1)
   {
      test_ = "SVecDVecSub - Addition assignment to right-hand side operand (1)";

      initialize();

      result_ =  db3_;
      result_ += sb3_ - db3_;
      db3_    += sb3_ - db3_;

      checkResult( db3_, result_ );
   }

   // Addition assignment to right-hand side operand (2)
   {
      test_ = "SVecDVecSub - Addition assignment to right-hand side operand (2)";

      initialize();

      result_ =  db3_;
      result_ += eval( sb3_ ) - db3_;
      db3_    += eval( sb3_ ) - db3_;

      checkResult( db3_, result_ );
   }

   // Addition assignment to right-hand side compound operand
   {
      test_ = "SVecDVecSub - Addition assignment to right-hand side compound operand";

      initialize();

      result_ =  db3_;
      result_ += sb3_ - ( dB3x3_ * db3_ );
      db3_    += sb3_ - ( dB3x3_ * db3_ );

      checkResult( db3_, result_ );
   }

   // Complex operation: a += ( 2*a ) - ( A * b );
   {
      test_ = "SVecDVecSub - Complex operation: a += ( 2*a ) - ( A * b );";

      initialize();

      result_ =  sb3_;
      result_ += ( 2*sb3_ ) - ( dA3x4_ * da4_ );
      sb3_    += ( 2*sb3_ ) - ( dA3x4_ * da4_ );

      checkResult( sb3_, result_ );
   }

   // Complex operation: a += ( A * b ) - ( 2*a );
   {
      test_ = "SVecDVecSub - Complex operation: a += ( A * b ) - ( 2*a );";

      initialize();

      result_ =  db3_;
      result_ += ( sA3x4_ * sa4_ ) - ( 2*db3_ );
      db3_    += ( sA3x4_ * sa4_ ) - ( 2*db3_ );

      checkResult( db3_, result_ );
   }


   //=====================================================================================
   // Subtraction with subtraction assignment
   //=====================================================================================

   // Subtraction assignment to left-hand side operand (1)
   {
      test_ = "SVecDVecSub - Subtraction assignment to left-hand side operand (1)";

      initialize();

      result_ =  sb3_;
      result_ -= sb3_ - db3_;
      sb3_    -= sb3_ - db3_;

      checkResult( sb3_, result_ );
   }

   // Subtraction assignment to left-hand side operand (2)
   {
      test_ = "SVecDVecSub - Subtraction assignment to left-hand side operand (2)";

      initialize();

      result_ =  sb3_;
      result_ -= sb3_ - eval( db3_ );
      sb3_    -= sb3_ - eval( db3_ );

      checkResult( sb3_, result_ );
   }

   // Subtraction assignment to left-hand side compound operand
   {
      test_ = "SVecDVecSub - Subtraction assignment to left-hand side compound operand";

      initialize();

      result_ =  sb3_;
      result_ -= ( sB3x3_ * sb3_ ) - db3_;
      sb3_    -= ( sB3x3_ * sb3_ ) - db3_;

      checkResult( sb3_, result_ );
   }

   // Subtraction assignment to right-hand side operand (1)
   {
      test_ = "SVecDVecSub - Subtraction assignment to right-hand side operand (1)";

      initialize();

      result_ =  db3_;
      result_ -= sb3_ - db3_;
      db3_    -= sb3_ - db3_;

      checkResult( db3_, result_ );
   }

   // Subtraction assignment to right-hand side operand (2)
   {
      test_ = "SVecDVecSub - Subtraction assignment to right-hand side operand (2)";

      initialize();

      result_ =  db3_;
      result_ -= eval( sb3_ ) - db3_;
      db3_    -= eval( sb3_ ) - db3_;

      checkResult( db3_, result_ );
   }

   // Subtraction assignment to right-hand side compound operand
   {
      test_ = "SVecDVecSub - Subtraction assignment to right-hand side compound operand";

      initialize();

      result_ =  db3_;
      result_ -= sb3_ - ( dB3x3_ * db3_ );
      db3_    -= sb3_ - ( dB3x3_ * db3_ );

      checkResult( db3_, result_ );
   }

   // Complex operation: a -= ( 2*a ) - ( A * b );
   {
      test_ = "SVecDVecSub - Complex operation: a -= ( 2*a ) - ( A * b );";

      initialize();

      result_ =  sb3_;
      result_ -= ( 2*sb3_ ) - ( dA3x4_ * da4_ );
      sb3_    -= ( 2*sb3_ ) - ( dA3x4_ * da4_ );

      checkResult( sb3_, result_ );
   }

   // Complex operation: a -= ( A * b ) - ( 2*a );
   {
      test_ = "SVecDVecSub - Complex operation: a -= ( A * b ) - ( 2*a );";

      initialize();

      result_ =  db3_;
      result_ -= ( sA3x4_ * sa4_ ) - ( 2*db3_ );
      db3_    -= ( sA3x4_ * sa4_ ) - ( 2*db3_ );

      checkResult( db3_, result_ );
   }


   //=====================================================================================
   // Subtraction with multiplication assignment
   //=====================================================================================

   // Multiplication assignment to left-hand side operand (1)
   {
      test_ = "SVecDVecSub - Multiplication assignment to left-hand side operand (1)";

      initialize();

      result_ =  sb3_;
      result_ *= sb3_ - db3_;
      sb3_    *= sb3_ - db3_;

      checkResult( sb3_, result_ );
   }

   // Multiplication assignment to left-hand side operand (2)
   {
      test_ = "SVecDVecSub - Multiplication assignment to left-hand side operand (2)";

      initialize();

      result_ =  sb3_;
      result_ *= sb3_ - eval( db3_ );
      sb3_    *= sb3_ - eval( db3_ );

      checkResult( sb3_, result_ );
   }

   // Multiplication assignment to left-hand side compound operand
   {
      test_ = "SVecDVecSub - Multiplication assignment to left-hand side compound operand";

      initialize();

      result_ =  sb3_;
      result_ *= ( sB3x3_ * sb3_ ) - db3_;
      sb3_    *= ( sB3x3_ * sb3_ ) - db3_;

      checkResult( sb3_, result_ );
   }

   // Multiplication assignment to right-hand side operand (1)
   {
      test_ = "SVecDVecSub - Multiplication assignment to right-hand side operand (1)";

      initialize();

      result_ =  db3_;
      result_ *= sb3_ - db3_;
      db3_    *= sb3_ - db3_;

      checkResult( db3_, result_ );
   }

   // Multiplication assignment to right-hand side operand (2)
   {
      test_ = "SVecDVecSub - Multiplication assignment to right-hand side operand (2)";

      initialize();

      result_ =  db3_;
      result_ *= eval( sb3_ ) - db3_;
      db3_    *= eval( sb3_ ) - db3_;

      checkResult( db3_, result_ );
   }

   // Multiplication assignment to right-hand side compound operand
   {
      test_ = "SVecDVecSub - Multiplication assignment to right-hand side compound operand";

      initialize();

      result_ =  db3_;
      result_ *= sb3_ - ( dB3x3_ * db3_ );
      db3_    *= sb3_ - ( dB3x3_ * db3_ );

      checkResult( db3_, result_ );
   }

   // Complex operation: a *= ( 2*a ) - ( A * b );
   {
      test_ = "SVecDVecSub - Complex operation: a *= ( 2*a ) - ( A * b );";

      initialize();

      result_ =  sb3_;
      result_ *= ( 2*sb3_ ) - ( dA3x4_ * da4_ );
      sb3_    *= ( 2*sb3_ ) - ( dA3x4_ * da4_ );

      checkResult( sb3_, result_ );
   }

   // Complex operation: a *= ( A * b ) - ( 2*a );
   {
      test_ = "SVecDVecSub - Complex operation: a *= ( A * b ) - ( 2*a );";

      initialize();

      result_ =  db3_;
      result_ *= ( sA3x4_ * sa4_ ) - ( 2*db3_ );
      db3_    *= ( sA3x4_ * sa4_ ) - ( 2*db3_ );

      checkResult( db3_, result_ );
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


   //=====================================================================================
   // Initialization of the dense vectors
   //=====================================================================================

   // Initializing the first dense column vector
   da4_.resize( 4UL, false );
   da4_[0] = -1;
   da4_[1] =  0;
   da4_[2] = -3;
   da4_[3] =  2;

   // Initializing the second dense column vector
   db3_.resize( 3UL, false );
   db3_[0] = 1;
   db3_[1] = 2;
   db3_[2] = 3;


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


   //=====================================================================================
   // Initialization of the dense matrices
   //=====================================================================================

   // Initializing the first row-major dense matrix
   dA3x4_(0,0) = -1;
   dA3x4_(0,1) =  0;
   dA3x4_(0,2) = -2;
   dA3x4_(0,3) =  0;
   dA3x4_(1,0) =  0;
   dA3x4_(1,1) =  2;
   dA3x4_(1,2) = -3;
   dA3x4_(1,3) =  1;
   dA3x4_(2,0) =  0;
   dA3x4_(2,1) =  1;
   dA3x4_(2,2) =  2;
   dA3x4_(2,3) =  2;

   // Initializing the second row-major dense matrix
   dB3x3_(0,0) =  0;
   dB3x3_(0,1) = -1;
   dB3x3_(0,2) =  0;
   dB3x3_(1,0) =  1;
   dB3x3_(1,1) = -2;
   dB3x3_(1,2) =  2;
   dB3x3_(2,0) =  0;
   dB3x3_(2,1) =  0;
   dB3x3_(2,2) = -3;
}
//*************************************************************************************************

} // namespace svecdvecsub

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
      RUN_SVECDVECSUB_ALIASING_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during aliasing test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
