//=================================================================================================
/*!
//  \file src/mathtest/operations/tdvecdmatmult/AliasingTest.cpp
//  \brief Source file for the dense vector/dense matrix multiplication aliasing test
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
#include <blazetest/mathtest/operations/tdvecdmatmult/AliasingTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace operations {

namespace tdvecdmatmult {

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
   : dA4x3_ ( 4UL, 3UL )
   , dB3x3_ ( 3UL, 3UL )
   , tdA4x3_( 4UL, 3UL )
   , tdB3x3_( 3UL, 3UL )
   , tda4_  ( 4UL )
   , tdb4_  ( 4UL )
   , tdc3_  ( 3UL )
   , tdd3_  ( 3UL )
   , tde3_  ( 3UL )
   , tsa4_  ( 4UL )
   , tsb3_  ( 3UL )
   , result_()
   , test_  ()
{
   testTDVecDMatMult ();
   testTDVecTDMatMult();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the dense vector/dense matrix multiplication.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs aliasing tests for the dense vector/dense matrix multiplication.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AliasingTest::testTDVecDMatMult()
{
   //=====================================================================================
   // Multiplication
   //=====================================================================================

   // Assignment to left-hand side operand
   {
      test_ = "TDVecDMatMult - Assignment to right-hand side vector operand";

      initialize();

      result_ = tda4_ * dA4x3_;
      tda4_   = tda4_ * dA4x3_;

      checkResult( tda4_, result_ );
   }

   // Assignment to first operand of left-hand side compound
   {
      test_ = "TDVecDMatMult - Assignment to first operand of left-hand side compound";

      initialize();

      result_ = tdb4_ * ( trans( tda4_ ) * tdc3_ );
      tda4_   = tdb4_ * ( trans( tda4_ ) * tdc3_ );

      checkResult( tda4_, result_ );
   }

   // Assignment to second operand of left-hand side compound
   {
      test_ = "TDVecDMatMult - Assignment to second operand of left-hand side compound";

      initialize();

      result_ = tdb4_ * ( trans( tda4_ ) * tdc3_ );
      tdc3_   = tdb4_ * ( trans( tda4_ ) * tdc3_ );

      checkResult( tdc3_, result_ );
   }

   // Assignment to first operand of right-hand side compound
   {
      test_ = "TDVecDMatMult - Assignment to first operand of right-hand side compound";

      initialize();

      result_ = ( tda4_ + tsa4_ ) * dA4x3_;
      tda4_   = ( tda4_ + tsa4_ ) * dA4x3_;

      checkResult( tda4_, result_ );
   }

   // Assignment to second operand of right-hand side compound
   {
      test_ = "TDVecDMatMult - Assignment to second operand of right-hand side compound";

      initialize();

      result_ = ( tda4_ + tsa4_ ) * dA4x3_;
      tsa4_   = ( tda4_ + tsa4_ ) * dA4x3_;

      checkResult( tsa4_, result_ );
   }


   //=====================================================================================
   // Multiplication with addition assignment
   //=====================================================================================

   // Addition assignment to left-hand side operand
   {
      test_ = "TDVecDMatMult - Addition assignment to right-hand side vector operand";

      initialize();

      result_ =  tdc3_;
      result_ += tdc3_ * dB3x3_;
      tdc3_   += tdc3_ * dB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Addition assignment to first operand of left-hand side compound
   {
      test_ = "TDVecDMatMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ += tde3_ * ( trans( tdc3_ ) * tdd3_ );
      tdc3_   += tde3_ * ( trans( tdc3_ ) * tdd3_ );

      checkResult( tdc3_, result_ );
   }

   // Addition assignment to second operand of left-hand side compound
   {
      test_ = "TDVecDMatMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tdd3_;
      result_ += tde3_ * ( trans( tdc3_ ) * tdd3_ );
      tdd3_   += tde3_ * ( trans( tdc3_ ) * tdd3_ );

      checkResult( tdd3_, result_ );
   }

   // Addition assignment to first operand of right-hand side compound
   {
      test_ = "TDVecDMatMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ += ( tdc3_ + tsb3_ ) * dB3x3_;
      tdc3_   += ( tdc3_ + tsb3_ ) * dB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Addition assignment to second operand of right-hand side compound
   {
      test_ = "TDVecDMatMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ += ( tdc3_ + tsb3_ ) * dB3x3_;
      tsb3_   += ( tdc3_ + tsb3_ ) * dB3x3_;

      checkResult( tsb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with subtraction assignment
   //=====================================================================================

   // Subtraction assignment to left-hand side operand
   {
      test_ = "TDVecDMatMult - Subtraction assignment to right-hand side vector operand";

      initialize();

      result_ =  tdc3_;
      result_ -= tdc3_ * dB3x3_;
      tdc3_   -= tdc3_ * dB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Subtraction assignment to first operand of left-hand side compound
   {
      test_ = "TDVecDMatMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ -= tde3_ * ( trans( tdc3_ ) * tdd3_ );
      tdc3_   -= tde3_ * ( trans( tdc3_ ) * tdd3_ );

      checkResult( tdc3_, result_ );
   }

   // Subtraction assignment to second operand of left-hand side compound
   {
      test_ = "TDVecDMatMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tdd3_;
      result_ -= tde3_ * ( trans( tdc3_ ) * tdd3_ );
      tdd3_   -= tde3_ * ( trans( tdc3_ ) * tdd3_ );

      checkResult( tdd3_, result_ );
   }

   // Subtraction assignment to first operand of right-hand side compound
   {
      test_ = "TDVecDMatMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ -= ( tdc3_ + tsb3_ ) * dB3x3_;
      tdc3_   -= ( tdc3_ + tsb3_ ) * dB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Subtraction assignment to second operand of right-hand side compound
   {
      test_ = "TDVecDMatMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ -= ( tdc3_ + tsb3_ ) * dB3x3_;
      tsb3_   -= ( tdc3_ + tsb3_ ) * dB3x3_;

      checkResult( tsb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with multiplication assignment
   //=====================================================================================

   // Multiplication assignment to left-hand side operand
   {
      test_ = "TDVecDMatMult - Multiplication assignment to right-hand side vector operand";

      initialize();

      result_ =  tdc3_;
      result_ *= tdc3_ * dB3x3_;
      tdc3_   *= tdc3_ * dB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Multiplication assignment to first operand of left-hand side compound
   {
      test_ = "TDVecDMatMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ *= tde3_ * ( trans( tdc3_ ) * tdd3_ );
      tdc3_   *= tde3_ * ( trans( tdc3_ ) * tdd3_ );

      checkResult( tdc3_, result_ );
   }

   // Multiplication assignment to second operand of left-hand side compound
   {
      test_ = "TDVecDMatMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tdd3_;
      result_ *= tde3_ * ( trans( tdc3_ ) * tdd3_ );
      tdd3_   *= tde3_ * ( trans( tdc3_ ) * tdd3_ );

      checkResult( tdd3_, result_ );
   }

   // Multiplication assignment to first operand of right-hand side compound
   {
      test_ = "TDVecDMatMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ *= ( tdc3_ + tsb3_ ) * dB3x3_;
      tdc3_   *= ( tdc3_ + tsb3_ ) * dB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Multiplication assignment to second operand of right-hand side compound
   {
      test_ = "TDVecDMatMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ *= ( tdc3_ + tsb3_ ) * dB3x3_;
      tsb3_   *= ( tdc3_ + tsb3_ ) * dB3x3_;

      checkResult( tsb3_, result_ );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the transpose dense vector/transpose dense matrix multiplication.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs aliasing tests for the transpose dense vector/transpose dense matrix
// multiplication. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AliasingTest::testTDVecTDMatMult()
{
   //=====================================================================================
   // Multiplication
   //=====================================================================================

   // Assignment to left-hand side operand
   {
      test_ = "TDVecTDMatMult - Assignment to right-hand side vector operand";

      initialize();

      result_ = tda4_ * tdA4x3_;
      tda4_   = tda4_ * tdA4x3_;

      checkResult( tda4_, result_ );
   }

   // Assignment to first operand of right-hand side compound
   {
      test_ = "TDVecTDMatMult - Assignment to first operand of right-hand side compound";

      initialize();

      result_ = ( tda4_ + tsa4_ ) * tdA4x3_;
      tda4_   = ( tda4_ + tsa4_ ) * tdA4x3_;

      checkResult( tda4_, result_ );
   }

   // Assignment to second operand of right-hand side compound
   {
      test_ = "TDVecTDMatMult - Assignment to second operand of right-hand side compound";

      initialize();

      result_ = ( tda4_ + tsa4_ ) * tdA4x3_;
      tsa4_   = ( tda4_ + tsa4_ ) * tdA4x3_;

      checkResult( tsa4_, result_ );
   }


   //=====================================================================================
   // Multiplication with addition assignment
   //=====================================================================================

   // Addition assignment to left-hand side operand
   {
      test_ = "TDVecTDMatMult - Addition assignment to right-hand side vector operand";

      initialize();

      result_ =  tdc3_;
      result_ += tdc3_ * tdB3x3_;
      tdc3_   += tdc3_ * tdB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Addition assignment to first operand of right-hand side compound
   {
      test_ = "TDVecTDMatMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ += ( tdc3_ + tsb3_ ) * tdB3x3_;
      tdc3_   += ( tdc3_ + tsb3_ ) * tdB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Addition assignment to second operand of right-hand side compound
   {
      test_ = "TDVecTDMatMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ += ( tdc3_ + tsb3_ ) * tdB3x3_;
      tsb3_   += ( tdc3_ + tsb3_ ) * tdB3x3_;

      checkResult( tsb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with subtraction assignment
   //=====================================================================================

   // Subtraction assignment to left-hand side operand
   {
      test_ = "TDVecTDMatMult - Subtraction assignment to right-hand side vector operand";

      initialize();

      result_ =  tdc3_;
      result_ -= tdc3_ * tdB3x3_;
      tdc3_   -= tdc3_ * tdB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Subtraction assignment to first operand of right-hand side compound
   {
      test_ = "TDVecTDMatMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ -= ( tdc3_ + tsb3_ ) * tdB3x3_;
      tdc3_   -= ( tdc3_ + tsb3_ ) * tdB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Subtraction assignment to second operand of right-hand side compound
   {
      test_ = "TDVecTDMatMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ -= ( tdc3_ + tsb3_ ) * tdB3x3_;
      tsb3_   -= ( tdc3_ + tsb3_ ) * tdB3x3_;

      checkResult( tsb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with multiplication assignment
   //=====================================================================================

   // Multiplication assignment to left-hand side operand
   {
      test_ = "TDVecTDMatMult - Multiplication assignment to right-hand side vector operand";

      initialize();

      result_ =  tdc3_;
      result_ *= tdc3_ * tdB3x3_;
      tdc3_   *= tdc3_ * tdB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Multiplication assignment to first operand of right-hand side compound
   {
      test_ = "TDVecTDMatMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tdc3_;
      result_ *= ( tdc3_ + tsb3_ ) * tdB3x3_;
      tdc3_   *= ( tdc3_ + tsb3_ ) * tdB3x3_;

      checkResult( tdc3_, result_ );
   }

   // Multiplication assignment to second operand of right-hand side compound
   {
      test_ = "TDVecTDMatMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tsb3_;
      result_ *= ( tdc3_ + tsb3_ ) * tdB3x3_;
      tsb3_   *= ( tdc3_ + tsb3_ ) * tdB3x3_;

      checkResult( tsb3_, result_ );
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
   // Initialization of the dense matrices
   //=====================================================================================

   // Initializing the first row-major dense matrix
   dA4x3_(0,0) = -1;
   dA4x3_(0,1) =  0;
   dA4x3_(0,2) = -2;
   dA4x3_(1,0) =  0;
   dA4x3_(1,1) =  2;
   dA4x3_(1,2) = -3;
   dA4x3_(2,0) =  0;
   dA4x3_(2,1) =  1;
   dA4x3_(2,2) =  2;
   dA4x3_(3,0) =  1;
   dA4x3_(3,1) =  0;
   dA4x3_(3,2) = -2;

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

   // Initializing the first column-major dense matrix
   tdA4x3_(0,0) = -1;
   tdA4x3_(0,1) =  0;
   tdA4x3_(0,2) = -2;
   tdA4x3_(1,0) =  0;
   tdA4x3_(1,1) =  2;
   tdA4x3_(1,2) = -3;
   tdA4x3_(2,0) =  0;
   tdA4x3_(2,1) =  1;
   tdA4x3_(2,2) =  2;
   tdA4x3_(3,0) =  1;
   tdA4x3_(3,1) =  0;
   tdA4x3_(3,2) = -2;

   // Initializing the second column-major dense matrix
   tdB3x3_(0,0) =  0;
   tdB3x3_(0,1) = -1;
   tdB3x3_(0,2) =  0;
   tdB3x3_(1,0) =  1;
   tdB3x3_(1,1) = -2;
   tdB3x3_(1,2) =  2;
   tdB3x3_(2,0) =  0;
   tdB3x3_(2,1) =  0;
   tdB3x3_(2,2) = -3;


   //=====================================================================================
   // Initialization of the dense vectors
   //=====================================================================================

   // Initializing the first dense row vector
   tda4_.resize( 4UL, false );
   tda4_[0] = -1;
   tda4_[1] =  0;
   tda4_[2] = -3;
   tda4_[3] =  2;

   // Initializing the second dense row vector
   tdb4_.resize( 4UL, false );
   tdb4_[0] =  0;
   tdb4_[1] =  1;
   tdb4_[2] =  2;
   tdb4_[3] = -1;

   // Initializing the third dense row vector
   tdc3_.resize( 3UL, false );
   tdc3_[0] = 1;
   tdc3_[1] = 2;
   tdc3_[2] = 3;

   // Initializing the fourth dense row vector
   tdd3_.resize( 3UL, false );
   tdd3_[0] = 0;
   tdd3_[1] = 2;
   tdd3_[2] = 1;

   // Initializing the fifth dense row vector
   tde3_.resize( 3UL, false );
   tde3_[0] = 0;
   tde3_[1] = 1;
   tde3_[2] = 3;


   //=====================================================================================
   // Initialization of the sparse vectors
   //=====================================================================================

   // Initializing the first sparse row vector
   tsa4_.resize( 4UL, false );
   tsa4_.reset();
   tsa4_[0] = -1;
   tsa4_[2] = -3;
   tsa4_[3] =  2;

   // Initializing the second sparse row vector
   tsb3_.resize( 3UL, false );
   tsb3_.reset();
   tsb3_[1] = 2;
   tsb3_[2] = 1;
}
//*************************************************************************************************

} // namespace tdvecdmatmult

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
      RUN_TDVECDMATMULT_ALIASING_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during aliasing test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
