//=================================================================================================
/*!
//  \file src/mathtest/operations/dmatdvecmult/AliasingTest.cpp
//  \brief Source file for the dense matrix/dense vector multiplication aliasing test
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
#include <blazetest/mathtest/operations/dmatdvecmult/AliasingTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace operations {

namespace dmatdvecmult {

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
   : dA3x4_ ( 3UL, 4UL )
   , dB3x3_ ( 3UL, 3UL )
   , tdA3x4_( 3UL, 4UL )
   , tdB3x3_( 3UL, 3UL )
   , da4_   ( 4UL )
   , db4_   ( 4UL )
   , dc3_   ( 3UL )
   , dd3_   ( 3UL )
   , de3_   ( 3UL )
   , sa4_   ( 4UL )
   , sb3_   ( 3UL )
   , result_()
   , test_  ()
{
   testDMatDVecMult ();
   testTDMatDVecMult();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the dense matrix/dense vector multiplication.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs aliasing tests for the dense matrix/dense vector multiplication.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AliasingTest::testDMatDVecMult()
{
   //=====================================================================================
   // Multiplication
   //=====================================================================================

   // Assignment to left-hand side operand
   {
      test_ = "DMatDVecMult - Assignment to right-hand side vector operand";

      initialize();

      result_ = dA3x4_ * da4_;
      da4_    = dA3x4_ * da4_;

      checkResult( da4_, result_ );
   }

   // Assignment to first operand of left-hand side compound
   {
      test_ = "DMatDVecMult - Assignment to first operand of left-hand side compound";

      initialize();

      result_ = ( dc3_ * trans( da4_ ) ) * db4_;
      dc3_    = ( dc3_ * trans( da4_ ) ) * db4_;

      checkResult( dc3_, result_ );
   }

   // Assignment to second operand of left-hand side compound
   {
      test_ = "DMatDVecMult - Assignment to second operand of left-hand side compound";

      initialize();

      result_ = ( dc3_ * trans( da4_ ) ) * db4_;
      da4_    = ( dc3_ * trans( da4_ ) ) * db4_;

      checkResult( da4_, result_ );
   }

   // Assignment to first operand of right-hand side compound
   {
      test_ = "DMatDVecMult - Assignment to first operand of right-hand side compound";

      initialize();

      result_ = dA3x4_ * ( da4_ + sa4_ );
      da4_    = dA3x4_ * ( da4_ + sa4_ );

      checkResult( da4_, result_ );
   }

   // Assignment to second operand of right-hand side compound
   {
      test_ = "DMatDVecMult - Assignment to second operand of right-hand side compound";

      initialize();

      result_ = dA3x4_ * ( da4_ + sa4_ );
      sa4_    = dA3x4_ * ( da4_ + sa4_ );

      checkResult( sa4_, result_ );
   }


   //=====================================================================================
   // Multiplication with addition assignment
   //=====================================================================================

   // Addition assignment to left-hand side operand
   {
      test_ = "DMatDVecMult - Addition assignment to right-hand side vector operand";

      initialize();

      result_ =  dc3_;
      result_ += dB3x3_ * dc3_;
      dc3_    += dB3x3_ * dc3_;

      checkResult( dc3_, result_ );
   }

   // Addition assignment to first operand of left-hand side compound
   {
      test_ = "DMatDVecMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  dc3_;
      result_ += ( dc3_ * trans( dd3_ ) ) * de3_;
      dc3_    += ( dc3_ * trans( dd3_ ) ) * de3_;

      checkResult( dc3_, result_ );
   }

   // Addition assignment to second operand of left-hand side compound
   {
      test_ = "DMatDVecMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  dd3_;
      result_ += ( dc3_ * trans( dd3_ ) ) * de3_;
      dd3_    += ( dc3_ * trans( dd3_ ) ) * de3_;

      checkResult( dd3_, result_ );
   }

   // Addition assignment to first operand of right-hand side compound
   {
      test_ = "DMatDVecMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  dc3_;
      result_ += dB3x3_ * ( dc3_ + sb3_ );
      dc3_    += dB3x3_ * ( dc3_ + sb3_ );

      checkResult( dc3_, result_ );
   }

   // Addition assignment to second operand of right-hand side compound
   {
      test_ = "DMatDVecMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  sb3_;
      result_ += dB3x3_ * ( dc3_ + sb3_ );
      sb3_    += dB3x3_ * ( dc3_ + sb3_ );

      checkResult( sb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with subtraction assignment
   //=====================================================================================

   // Subtraction assignment to left-hand side operand
   {
      test_ = "DMatDVecMult - Subtraction assignment to right-hand side vector operand";

      initialize();

      result_ =  dc3_;
      result_ -= dB3x3_ * dc3_;
      dc3_    -= dB3x3_ * dc3_;

      checkResult( dc3_, result_ );
   }

   // Subtraction assignment to first operand of left-hand side compound
   {
      test_ = "DMatDVecMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  dc3_;
      result_ -= ( dc3_ * trans( dd3_ ) ) * de3_;
      dc3_    -= ( dc3_ * trans( dd3_ ) ) * de3_;

      checkResult( dc3_, result_ );
   }

   // Subtraction assignment to second operand of left-hand side compound
   {
      test_ = "DMatDVecMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  dd3_;
      result_ -= ( dc3_ * trans( dd3_ ) ) * de3_;
      dd3_    -= ( dc3_ * trans( dd3_ ) ) * de3_;

      checkResult( dd3_, result_ );
   }

   // Subtraction assignment to first operand of right-hand side compound
   {
      test_ = "DMatDVecMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  dc3_;
      result_ -= dB3x3_ * ( dc3_ + sb3_ );
      dc3_    -= dB3x3_ * ( dc3_ + sb3_ );

      checkResult( dc3_, result_ );
   }

   // Subtraction assignment to second operand of right-hand side compound
   {
      test_ = "DMatDVecMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  sb3_;
      result_ -= dB3x3_ * ( dc3_ + sb3_ );
      sb3_    -= dB3x3_ * ( dc3_ + sb3_ );

      checkResult( sb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with multiplication assignment
   //=====================================================================================

   // Multiplication assignment to left-hand side operand
   {
      test_ = "DMatDVecMult - Multiplication assignment to right-hand side vector operand";

      initialize();

      result_ =  dc3_;
      result_ *= dB3x3_ * dc3_;
      dc3_    *= dB3x3_ * dc3_;

      checkResult( dc3_, result_ );
   }

   // Multiplication assignment to first operand of left-hand side compound
   {
      test_ = "DMatDVecMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  dc3_;
      result_ *= ( dc3_ * trans( dd3_ ) ) * de3_;
      dc3_    *= ( dc3_ * trans( dd3_ ) ) * de3_;

      checkResult( dc3_, result_ );
   }

   // Multiplication assignment to second operand of left-hand side compound
   {
      test_ = "DMatDVecMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  dd3_;
      result_ *= ( dc3_ * trans( dd3_ ) ) * de3_;
      dd3_    *= ( dc3_ * trans( dd3_ ) ) * de3_;

      checkResult( dd3_, result_ );
   }

   // Multiplication assignment to first operand of right-hand side compound
   {
      test_ = "DMatDVecMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  dc3_;
      result_ *= dB3x3_ * ( dc3_ + sb3_ );
      dc3_    *= dB3x3_ * ( dc3_ + sb3_ );

      checkResult( dc3_, result_ );
   }

   // Multiplication assignment to second operand of right-hand side compound
   {
      test_ = "DMatDVecMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  sb3_;
      result_ *= dB3x3_ * ( dc3_ + sb3_ );
      sb3_    *= dB3x3_ * ( dc3_ + sb3_ );

      checkResult( sb3_, result_ );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the transpose dense matrix/dense vector multiplication.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs aliasing tests for the transpose dense matrix/dense vector
// multiplication. In case an error is detected, a \a std::runtime_error exception is
// thrown.
*/
void AliasingTest::testTDMatDVecMult()
{
   //=====================================================================================
   // Multiplication
   //=====================================================================================

   // Assignment to left-hand side operand
   {
      test_ = "TDMatDVecMult - Assignment to right-hand side vector operand";

      initialize();

      result_ = tdA3x4_ * da4_;
      da4_    = tdA3x4_ * da4_;

      checkResult( da4_, result_ );
   }

   // Assignment to first operand of right-hand side compound
   {
      test_ = "TDMatDVecMult - Assignment to first operand of right-hand side compound";

      initialize();

      result_ = tdA3x4_ * ( da4_ + sa4_ );
      da4_    = tdA3x4_ * ( da4_ + sa4_ );

      checkResult( da4_, result_ );
   }

   // Assignment to second operand of right-hand side compound
   {
      test_ = "TDMatDVecMult - Assignment to second operand of right-hand side compound";

      initialize();

      result_ = tdA3x4_ * ( da4_ + sa4_ );
      sa4_    = tdA3x4_ * ( da4_ + sa4_ );

      checkResult( sa4_, result_ );
   }


   //=====================================================================================
   // Multiplication with addition assignment
   //=====================================================================================

   // Addition assignment to left-hand side operand
   {
      test_ = "TDMatDVecMult - Addition assignment to right-hand side vector operand";

      initialize();

      result_ =  dc3_;
      result_ += tdB3x3_ * dc3_;
      dc3_    += tdB3x3_ * dc3_;

      checkResult( dc3_, result_ );
   }

   // Addition assignment to first operand of right-hand side compound
   {
      test_ = "TDMatDVecMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  dc3_;
      result_ += tdB3x3_ * ( dc3_ + sb3_ );
      dc3_    += tdB3x3_ * ( dc3_ + sb3_ );

      checkResult( dc3_, result_ );
   }

   // Addition assignment to second operand of right-hand side compound
   {
      test_ = "TDMatDVecMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  sb3_;
      result_ += tdB3x3_ * ( dc3_ + sb3_ );
      sb3_    += tdB3x3_ * ( dc3_ + sb3_ );

      checkResult( sb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with subtraction assignment
   //=====================================================================================

   // Subtraction assignment to left-hand side operand
   {
      test_ = "TDMatDVecMult - Subtraction assignment to right-hand side vector operand";

      initialize();

      result_ =  dc3_;
      result_ -= tdB3x3_ * dc3_;
      dc3_    -= tdB3x3_ * dc3_;

      checkResult( dc3_, result_ );
   }

   // Subtraction assignment to first operand of right-hand side compound
   {
      test_ = "TDMatDVecMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  dc3_;
      result_ -= tdB3x3_ * ( dc3_ + sb3_ );
      dc3_    -= tdB3x3_ * ( dc3_ + sb3_ );

      checkResult( dc3_, result_ );
   }

   // Subtraction assignment to second operand of right-hand side compound
   {
      test_ = "TDMatDVecMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  sb3_;
      result_ -= tdB3x3_ * ( dc3_ + sb3_ );
      sb3_    -= tdB3x3_ * ( dc3_ + sb3_ );

      checkResult( sb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with multiplication assignment
   //=====================================================================================

   // Multiplication assignment to left-hand side operand
   {
      test_ = "TDMatDVecMult - Multiplication assignment to right-hand side vector operand";

      initialize();

      result_ =  dc3_;
      result_ *= tdB3x3_ * dc3_;
      dc3_    *= tdB3x3_ * dc3_;

      checkResult( dc3_, result_ );
   }

   // Multiplication assignment to first operand of right-hand side compound
   {
      test_ = "TDMatDVecMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  dc3_;
      result_ *= tdB3x3_ * ( dc3_ + sb3_ );
      dc3_    *= tdB3x3_ * ( dc3_ + sb3_ );

      checkResult( dc3_, result_ );
   }

   // Multiplication assignment to second operand of right-hand side compound
   {
      test_ = "TDMatDVecMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  sb3_;
      result_ *= tdB3x3_ * ( dc3_ + sb3_ );
      sb3_    *= tdB3x3_ * ( dc3_ + sb3_ );

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

   // Initializing the first column-major dense matrix
   tdA3x4_(0,0) = -1;
   tdA3x4_(0,1) =  0;
   tdA3x4_(0,2) = -2;
   tdA3x4_(0,3) =  0;
   tdA3x4_(1,0) =  0;
   tdA3x4_(1,1) =  2;
   tdA3x4_(1,2) = -3;
   tdA3x4_(1,3) =  1;
   tdA3x4_(2,0) =  0;
   tdA3x4_(2,1) =  1;
   tdA3x4_(2,2) =  2;
   tdA3x4_(2,3) =  2;

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

   // Initializing the first dense column vector
   da4_.resize( 4UL, false );
   da4_[0] = -1;
   da4_[1] =  0;
   da4_[2] = -3;
   da4_[3] =  2;

   // Initializing the second dense column vector
   db4_.resize( 4UL, false );
   db4_[0] =  0;
   db4_[1] =  1;
   db4_[2] =  2;
   db4_[3] = -1;

   // Initializing the third dense column vector
   dc3_.resize( 3UL, false );
   dc3_[0] = 1;
   dc3_[1] = 2;
   dc3_[2] = 3;

   // Initializing the fourth dense column vector
   dd3_.resize( 3UL, false );
   dd3_[0] = 0;
   dd3_[1] = 2;
   dd3_[2] = 1;

   // Initializing the fifth dense column vector
   de3_.resize( 3UL, false );
   de3_[0] = 0;
   de3_[1] = 1;
   de3_[2] = 3;


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
   sb3_[1] = 2;
   sb3_[2] = 1;
}
//*************************************************************************************************

} // namespace dmatdvecmult

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
      RUN_DMATDVECMULT_ALIASING_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during aliasing test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
