//=================================================================================================
/*!
//  \file src/mathtest/operations/tsvecsmatmult/AliasingTest.cpp
//  \brief Source file for the sparse vector/sparse matrix multiplication aliasing test
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
#include <blazetest/mathtest/operations/tsvecsmatmult/AliasingTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace operations {

namespace tsvecsmatmult {

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
   : sA4x3_ ( 4UL, 3UL )
   , sB3x3_ ( 3UL, 3UL )
   , tsA4x3_( 4UL, 3UL )
   , tsB3x3_( 3UL, 3UL )
   , tsa4_  ( 4UL )
   , tsb4_  ( 4UL )
   , tsc3_  ( 3UL )
   , tsd3_  ( 3UL )
   , tda4_  ( 4UL )
   , tdb3_  ( 3UL )
   , result_()
   , test_  ()
{
   testTSVecSMatMult ();
   testTSVecTSMatMult();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the sparse vector/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs aliasing tests for the sparse vector/sparse matrix multiplication.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AliasingTest::testTSVecSMatMult()
{
   //=====================================================================================
   // Multiplication
   //=====================================================================================

   // Assignment to left-hand side operand
   {
      test_ = "TSVecSMatMult - Assignment to right-hand side vector operand";

      initialize();

      result_ = tsa4_ * sA4x3_;
      tsa4_   = tsa4_ * sA4x3_;

      checkResult( tsa4_, result_ );
   }

   // Assignment to first operand of left-hand side compound
   {
      test_ = "TSVecSMatMult - Assignment to first operand of left-hand side compound";

      initialize();

      result_ = tsb4_ * ( trans( tsa4_ ) * tdb3_ );
      tsa4_   = tsb4_ * ( trans( tsa4_ ) * tdb3_ );

      checkResult( tsa4_, result_ );
   }

   // Assignment to second operand of left-hand side compound
   {
      test_ = "TSVecSMatMult - Assignment to second operand of left-hand side compound";

      initialize();

      result_ = tsb4_ * ( trans( tsa4_ ) * tdb3_ );
      tdb3_   = tsb4_ * ( trans( tsa4_ ) * tdb3_ );

      checkResult( tdb3_, result_ );
   }

   // Assignment to first operand of right-hand side compound
   {
      test_ = "TSVecSMatMult - Assignment to first operand of right-hand side compound";

      initialize();

      result_ = ( tsa4_ * tda4_ ) * sA4x3_;
      tsa4_   = ( tsa4_ * tda4_ ) * sA4x3_;

      checkResult( tsa4_, result_ );
   }

   // Assignment to second operand of right-hand side compound
   {
      test_ = "TSVecSMatMult - Assignment to second operand of right-hand side compound";

      initialize();

      result_ = ( tsa4_ * tda4_ ) * sA4x3_;
      tda4_   = ( tsa4_ * tda4_ ) * sA4x3_;

      checkResult( tda4_, result_ );
   }


   //=====================================================================================
   // Multiplication with addition assignment
   //=====================================================================================

   // Addition assignment to left-hand side operand
   {
      test_ = "TSVecSMatMult - Addition assignment to right-hand side vector operand";

      initialize();

      result_ =  tsc3_;
      result_ += tsc3_ * sB3x3_;
      tsc3_   += tsc3_ * sB3x3_;

      checkResult( tsc3_, result_ );
   }

   // Addition assignment to first operand of left-hand side compound
   {
      test_ = "TSVecSMatMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tsc3_;
      result_ += tsd3_ * ( trans( tsc3_ ) * tdb3_ );
      tsc3_   += tsd3_ * ( trans( tsc3_ ) * tdb3_ );

      checkResult( tsc3_, result_ );
   }

   // Addition assignment to second operand of left-hand side compound
   {
      test_ = "TSVecSMatMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tdb3_;
      result_ += tsd3_ * ( trans( tdb3_ ) * tsc3_ );
      tdb3_   += tsd3_ * ( trans( tdb3_ ) * tsc3_ );

      checkResult( tdb3_, result_ );
   }

   // Addition assignment to first operand of right-hand side compound
   {
      test_ = "TSVecSMatMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tsc3_;
      result_ += ( tsc3_ * tdb3_ ) * sB3x3_;
      tsc3_   += ( tsc3_ * tdb3_ ) * sB3x3_;

      checkResult( tsc3_, result_ );
   }

   // Addition assignment to second operand of right-hand side compound
   {
      test_ = "TSVecSMatMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tdb3_;
      result_ += ( tsc3_ * tdb3_ ) * sB3x3_;
      tdb3_   += ( tsc3_ * tdb3_ ) * sB3x3_;

      checkResult( tdb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with subtraction assignment
   //=====================================================================================

   // Subtraction assignment to left-hand side operand
   {
      test_ = "TSVecSMatMult - Subtraction assignment to right-hand side vector operand";

      initialize();

      result_ =  tsc3_;
      result_ -= tsc3_ * sB3x3_;
      tsc3_   -= tsc3_ * sB3x3_;

      checkResult( tsc3_, result_ );
   }

   // Subtraction assignment to first operand of left-hand side compound
   {
      test_ = "TSVecSMatMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tsc3_;
      result_ -= tsd3_ * ( trans( tsc3_ ) * tdb3_ );
      tsc3_   -= tsd3_ * ( trans( tsc3_ ) * tdb3_ );

      checkResult( tsc3_, result_ );
   }

   // Subtraction assignment to second operand of left-hand side compound
   {
      test_ = "TSVecSMatMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tdb3_;
      result_ -= tsd3_ * ( trans( tdb3_ ) * tsc3_ );
      tdb3_   -= tsd3_ * ( trans( tdb3_ ) * tsc3_ );

      checkResult( tdb3_, result_ );
   }

   // Subtraction assignment to first operand of right-hand side compound
   {
      test_ = "TSVecSMatMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tsc3_;
      result_ -= ( tsc3_ * tdb3_ ) * sB3x3_;
      tsc3_   -= ( tsc3_ * tdb3_ ) * sB3x3_;

      checkResult( tsc3_, result_ );
   }

   // Subtraction assignment to second operand of right-hand side compound
   {
      test_ = "TSVecSMatMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tdb3_;
      result_ -= ( tsc3_ * tdb3_ ) * sB3x3_;
      tdb3_   -= ( tsc3_ * tdb3_ ) * sB3x3_;

      checkResult( tdb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with multiplication assignment
   //=====================================================================================

   // Multiplication assignment to left-hand side operand
   {
      test_ = "TSVecSMatMult - Multiplication assignment to right-hand side vector operand";

      initialize();

      result_ =  tsc3_;
      result_ *= tsc3_ * sB3x3_;
      tsc3_   *= tsc3_ * sB3x3_;

      checkResult( tsc3_, result_ );
   }

   // Multiplication assignment to first operand of left-hand side compound
   {
      test_ = "TSVecSMatMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tsc3_;
      result_ *= tsd3_ * ( trans( tsc3_ ) * tdb3_ );
      tsc3_   *= tsd3_ * ( trans( tsc3_ ) * tdb3_ );

      checkResult( tsc3_, result_ );
   }

   // Multiplication assignment to second operand of left-hand side compound
   {
      test_ = "TSVecSMatMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tdb3_;
      result_ *= tsd3_ * ( trans( tdb3_ ) * tsc3_ );
      tdb3_   *= tsd3_ * ( trans( tdb3_ ) * tsc3_ );

      checkResult( tdb3_, result_ );
   }

   // Multiplication assignment to first operand of right-hand side compound
   {
      test_ = "TSVecSMatMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tsc3_;
      result_ *= ( tsc3_ * tdb3_ ) * sB3x3_;
      tsc3_   *= ( tsc3_ * tdb3_ ) * sB3x3_;

      checkResult( tsc3_, result_ );
   }

   // Multiplication assignment to second operand of right-hand side compound
   {
      test_ = "TSVecSMatMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tdb3_;
      result_ *= ( tsc3_ * tdb3_ ) * sB3x3_;
      tdb3_   *= ( tsc3_ * tdb3_ ) * sB3x3_;

      checkResult( tdb3_, result_ );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the transpose sparse vector/transpose sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs aliasing tests for the transpose sparse vector/transpose sparse matrix
// multiplication. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void AliasingTest::testTSVecTSMatMult()
{
   //=====================================================================================
   // Multiplication
   //=====================================================================================

   // Assignment to left-hand side operand
   {
      test_ = "TSVecTSMatMult - Assignment to right-hand side vector operand";

      initialize();

      result_ = tsa4_ * tsA4x3_;
      tsa4_   = tsa4_ * tsA4x3_;

      checkResult( tsa4_, result_ );
   }

   // Assignment to first operand of right-hand side compound
   {
      test_ = "TSVecTSMatMult - Assignment to first operand of right-hand side compound";

      initialize();

      result_ = ( tsa4_ * tda4_ ) * tsA4x3_;
      tsa4_   = ( tsa4_ * tda4_ ) * tsA4x3_;

      checkResult( tsa4_, result_ );
   }

   // Assignment to second operand of right-hand side compound
   {
      test_ = "TSVecTSMatMult - Assignment to second operand of right-hand side compound";

      initialize();

      result_ = ( tsa4_ * tda4_ ) * tsA4x3_;
      tda4_   = ( tsa4_ * tda4_ ) * tsA4x3_;

      checkResult( tda4_, result_ );
   }


   //=====================================================================================
   // Multiplication with addition assignment
   //=====================================================================================

   // Addition assignment to left-hand side operand
   {
      test_ = "TSVecTSMatMult - Addition assignment to right-hand side vector operand";

      initialize();

      result_ =  tsc3_;
      result_ += tsc3_ * tsB3x3_;
      tsc3_   += tsc3_ * tsB3x3_;

      checkResult( tsc3_, result_ );
   }

   // Addition assignment to first operand of right-hand side compound
   {
      test_ = "TSVecTSMatMult - Addition assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tsc3_;
      result_ += ( tsc3_ * tdb3_ ) * tsB3x3_;
      tsc3_   += ( tsc3_ * tdb3_ ) * tsB3x3_;

      checkResult( tsc3_, result_ );
   }

   // Addition assignment to second operand of right-hand side compound
   {
      test_ = "TSVecTSMatMult - Addition assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tdb3_;
      result_ += ( tsc3_ * tdb3_ ) * tsB3x3_;
      tdb3_   += ( tsc3_ * tdb3_ ) * tsB3x3_;

      checkResult( tdb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with subtraction assignment
   //=====================================================================================

   // Subtraction assignment to left-hand side operand
   {
      test_ = "TSVecTSMatMult - Subtraction assignment to right-hand side vector operand";

      initialize();

      result_ =  tsc3_;
      result_ -= tsc3_ * tsB3x3_;
      tsc3_   -= tsc3_ * tsB3x3_;

      checkResult( tsc3_, result_ );
   }

   // Subtraction assignment to first operand of right-hand side compound
   {
      test_ = "TSVecTSMatMult - Subtraction assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tsc3_;
      result_ -= ( tsc3_ * tdb3_ ) * tsB3x3_;
      tsc3_   -= ( tsc3_ * tdb3_ ) * tsB3x3_;

      checkResult( tsc3_, result_ );
   }

   // Subtraction assignment to second operand of right-hand side compound
   {
      test_ = "TSVecTSMatMult - Subtraction assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tdb3_;
      result_ -= ( tsc3_ * tdb3_ ) * tsB3x3_;
      tdb3_   -= ( tsc3_ * tdb3_ ) * tsB3x3_;

      checkResult( tdb3_, result_ );
   }


   //=====================================================================================
   // Multiplication with multiplication assignment
   //=====================================================================================

   // Multiplication assignment to left-hand side operand
   {
      test_ = "TSVecTSMatMult - Multiplication assignment to right-hand side vector operand";

      initialize();

      result_ =  tsc3_;
      result_ *= tsc3_ * tsB3x3_;
      tsc3_   *= tsc3_ * tsB3x3_;

      checkResult( tsc3_, result_ );
   }

   // Multiplication assignment to first operand of right-hand side compound
   {
      test_ = "TSVecTSMatMult - Multiplication assignment to first operand of left-hand side compound";

      initialize();

      result_ =  tsc3_;
      result_ *= ( tsc3_ * tdb3_ ) * tsB3x3_;
      tsc3_   *= ( tsc3_ * tdb3_ ) * tsB3x3_;

      checkResult( tsc3_, result_ );
   }

   // Multiplication assignment to second operand of right-hand side compound
   {
      test_ = "TSVecTSMatMult - Multiplication assignment to second operand of left-hand side compound";

      initialize();

      result_ =  tdb3_;
      result_ *= ( tsc3_ * tdb3_ ) * tsB3x3_;
      tdb3_   *= ( tsc3_ * tdb3_ ) * tsB3x3_;

      checkResult( tdb3_, result_ );
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
   // Initialization of the sparse matrices
   //=====================================================================================

   // Initializing the first row-major sparse matrix
   sA4x3_(0,0) = -1;
   sA4x3_(0,2) = -2;
   sA4x3_(1,1) =  2;
   sA4x3_(1,2) = -3;
   sA4x3_(2,1) =  1;
   sA4x3_(2,2) =  2;
   sA4x3_(3,0) =  1;
   sA4x3_(3,2) = -2;

   // Initializing the second row-major sparse matrix
   sB3x3_(0,0) = -1;
   sB3x3_(1,0) =  1;
   sB3x3_(1,1) = -2;
   sB3x3_(1,2) =  2;
   sB3x3_(2,2) = -3;

   // Initializing the first column-major sparse matrix
   tsA4x3_(0,0) = -1;
   tsA4x3_(0,2) = -2;
   tsA4x3_(1,1) =  2;
   tsA4x3_(1,2) = -3;
   tsA4x3_(2,1) =  1;
   tsA4x3_(2,2) =  2;
   tsA4x3_(3,0) =  1;
   tsA4x3_(3,2) = -2;

   // Initializing the second column-major sparse matrix
   tsB3x3_(0,0) = -1;
   tsB3x3_(1,0) =  1;
   tsB3x3_(1,1) = -2;
   tsB3x3_(1,2) =  2;
   tsB3x3_(2,2) = -3;


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
   tsb4_.resize( 4UL, false );
   tsb4_.reset();
   tsb4_[1] =  1;
   tsb4_[2] =  2;
   tsb4_[3] = -1;

   // Initializing the third sparse row vector
   tsc3_.resize( 3UL, false );
   tsc3_.reset();
   tsc3_[0] = 1;
   tsc3_[1] = 2;
   tsc3_[2] = 3;

   // Initializing the fourth sparse row vector
   tsd3_.resize( 3UL, false );
   tsd3_.reset();
   tsd3_[1] = 2;
   tsd3_[2] = 1;


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
   tdb3_.resize( 3UL, false );
   tdb3_[0] = 0;
   tdb3_[1] = 2;
   tdb3_[2] = 1;
}
//*************************************************************************************************

} // namespace tsvecsmatmult

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
      RUN_TSVECSMATMULT_ALIASING_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during aliasing test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
