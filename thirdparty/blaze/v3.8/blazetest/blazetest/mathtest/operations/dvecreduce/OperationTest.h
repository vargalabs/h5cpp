//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/dvecreduce/OperationTest.h
//  \brief Header file for the dense vector total reduction operation test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_DVECREDUCE_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_DVECREDUCE_OPERATIONTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/Aliases.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/functors/Add.h>
#include <blaze/math/traits/ReduceTrait.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/Random.h>
#include <blazetest/system/MathTest.h>
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/IsEqual.h>


namespace blazetest {

namespace mathtest {

namespace operations {

namespace dvecreduce {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the dense vector reduction operation test.
//
// This class template represents one particular test of a reduction operation on a vector of a
// particular type. The template argument \a VT represents the type of the vector operand.
*/
template< typename VT >  // Type of the dense vector
class OperationTest
{
 private:
   //**Type definitions****************************************************************************
   using ET = blaze::ElementType_t<VT>;  //!< Element type.

   using TVT = blaze::TransposeType_t<VT>;  //!< Transpose vector type.

   using RE = blaze::ReduceTrait_t<VT,blaze::Add>;  //!< Result type of the reduction operation.

   using RT = blaze::CompressedVector<ET,false>;  //!< Reference type.
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename OP >
   explicit OperationTest( const Creator<VT>& creator, OP op );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

 private:
   //**Test functions******************************************************************************
   /*!\name Test functions */
   //@{
                           void testInitialStatus ();
                           void testAssignment    ();
   template< typename OP > void testBasicOperation( OP op );
   //@}
   //**********************************************************************************************

   //**Error detection functions*******************************************************************
   /*!\name Error detection functions */
   //@{
   void checkResult();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initResults();
   void convertException( const std::exception& ex );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VT vec_;     //!< The dense vector operand.
   RE res_;     //!< The result of the reduction operation.
   RT refvec_;  //!< The reference vector.
   RE refres_;  //!< The reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TVT );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( RT  );

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( VT  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( TVT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( RT  );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET, blaze::ElementType_t<TVT> );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the dense vector reduction operation test.
//
// \param creator The creator for dense vector operand.
// \param op The reduction operation.
// \exception std::runtime_error Operation error detected.
*/
template< typename VT >  // Type of the dense vector
template< typename OP >  // Type of the reduction operation
OperationTest<VT>::OperationTest( const Creator<VT>& creator, OP op )
   : vec_( creator( NoZeros() ) )  // The dense vector operand
   , res_()                        // The result of the reduction operation
   , refvec_( vec_ )               // The reference vector
   , refres_()                     // The reference result
   , test_()                       // Label of the currently performed test
   , error_()                      // Description of the current error type
{
   testInitialStatus();
   testAssignment();
   testBasicOperation( op );
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Tests on the initial status of the vector operand.
//
// \return void
// \exception std::runtime_error Initialization error detected.
//
// This function runs tests on the initial status of the vector operand. In case any initialization
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT >  // Type of the dense vector
void OperationTest<VT>::testInitialStatus()
{
   // Checking the size of the dense operand
   if( vec_.size() != refvec_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of dense operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Detected size = " << vec_.size() << "\n"
          << "   Expected size = " << refvec_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the dense operand
   if( !isEqual( vec_, refvec_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of dense operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Current initialization:\n" << vec_ << "\n"
          << "   Expected initialization:\n" << refvec_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the vector assignment.
//
// \return void
// \exception std::runtime_error Assignment error detected.
//
// This function tests the vector assignment. In case any error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT >  // Type of the dense vector
void OperationTest<VT>::testAssignment()
{
   try {
      vec_ = refvec_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the given vector\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( vec_, refvec_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of dense operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Current initialization:\n" << vec_ << "\n"
          << "   Expected initialization:\n" << refvec_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain dense vector reduction operation.
//
// \param op The reduction operation.
// \return void
// \exception std::runtime_error Reduction error detected.
//
// This function tests the plain reduction operation with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// reduction or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT >  // Type of the left-hand side dense vector
template< typename OP >  // Type of the reduction operation
void OperationTest<VT>::testBasicOperation( OP op )
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Reduction operation
      //=====================================================================================

      // Reduction operation with the given vector
      {
         test_  = "Reduction operation with the given vector";
         error_ = "Failed reduction operation";

         try {
            initResults();
            res_    = reduce( vec_, op );
            refres_ = reduce( refvec_, op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }

      // Reduction operation with evaluated vector
      {
         test_  = "Reduction operation with evaluated vector";
         error_ = "Failed reduction operation";

         try {
            initResults();
            res_    = reduce( eval( vec_ ), op );
            refres_ = reduce( eval( refvec_ ), op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }


      //=====================================================================================
      // Reduction operation with addition assignment
      //=====================================================================================

      // Reduction operation with addition assignment with the given vectors
      {
         test_  = "Reduction operation with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            res_    += reduce( vec_, op );
            refres_ += reduce( refvec_, op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }

      // Reduction operation with addition assignment with evaluated vector
      {
         test_  = "Reduction operation with addition assignment with evaluated vector";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            res_    += reduce( eval( vec_ ), op );
            refres_ += reduce( eval( refvec_ ), op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }


      //=====================================================================================
      // Reduction operation with subtraction assignment
      //=====================================================================================

      // Reduction operation with subtraction assignment with the given vectors
      {
         test_  = "Reduction operation with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            res_    -= reduce( vec_, op );
            refres_ -= reduce( refvec_, op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }

      // Reduction operation with subtraction assignment with evaluated vector
      {
         test_  = "Reduction operation with subtraction assignment with evaluated vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            res_    -= reduce( eval( vec_ ), op );
            refres_ -= reduce( eval( refvec_ ), op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }


      //=====================================================================================
      // Reduction operation with multiplication assignment
      //=====================================================================================

      // Reduction operation with multiplication assignment with the given vectors
      {
         test_  = "Reduction operation with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            res_    *= reduce( vec_, op );
            refres_ *= reduce( refvec_, op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }

      // Reduction operation with multiplication assignment with evaluated vector
      {
         test_  = "Reduction operation with multiplication assignment with evaluated vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            res_    *= reduce( eval( vec_ ), op );
            refres_ *= reduce( eval( refvec_ ), op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }
   }
#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  ERROR DETECTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checking and comparing the computed results.
//
// \return void
// \exception std::runtime_error Incorrect result detected.
//
// This function is called after each test case to check and compare the computed results.
*/
template< typename VT >  // Type of the dense vector
void OperationTest<VT>::checkResult()
{
   if( !isEqual( res_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Result:\n" << res_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initializing the results.
//
// \return void
//
// This function is called before each test case to initialize the results to random values.
*/
template< typename VT >  // Type of the dense vector
void OperationTest<VT>::initResults()
{
   using blaze::randomize;

   const blaze::UnderlyingBuiltin_t<RE> min( randmin );
   const blaze::UnderlyingBuiltin_t<RE> max( randmax );

   randomize( res_, min, max );

   refres_ = res_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Convert the given exception into a \a std::runtime_error exception.
//
// \param ex The \a std::exception to be extended.
// \return void
// \exception std::runtime_error The converted exception.
//
// This function converts the given exception to a \a std::runtime_error exception. Additionally,
// the function extends the given exception message by all available information for the failed
// test.
*/
template< typename VT >  // Type of the dense vector
void OperationTest<VT>::convertException( const std::exception& ex )
{
   std::ostringstream oss;
   oss << " Test : " << test_ << "\n"
       << " Error: " << error_ << "\n"
       << " Details:\n"
       << "   Random seed = " << blaze::getSeed() << "\n"
       << "   Dense vector type:\n"
       << "     " << typeid( VT ).name() << "\n"
       << "   Error message: " << ex.what() << "\n";
   throw std::runtime_error( oss.str() );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the reduction operation for a specific vector type.
//
// \param creator The creator for the dense vector.
// \return void
*/
template< typename VT >  // Type of the dense vector
void runTest( const Creator<VT>& creator )
{
#if BLAZETEST_MATHTEST_TEST_ADDITION
   if( BLAZETEST_MATHTEST_TEST_ADDITION > 1 )
   {
      class Sum : public blaze::Add
      {};

      for( size_t rep=0UL; rep<BLAZETEST_REPETITIONS; ++rep ) {
         OperationTest<VT>( creator, []( const auto& a, const auto& b ){ return a + b; } );
         OperationTest<VT>( creator, blaze::Add() );
         OperationTest<VT>( creator, Sum() );
      }
   }
#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  MACROS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the definition of a dense vector reduction operation test case.
*/
#define DEFINE_DVECREDUCE_OPERATION_TEST( VT ) \
   extern template class blazetest::mathtest::operations::dvecreduce::OperationTest<VT>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a dense vector reduction operation test case.
*/
#define RUN_DVECREDUCE_OPERATION_TEST( C ) \
   blazetest::mathtest::operations::dvecreduce::runTest( C )
/*! \endcond */
//*************************************************************************************************

} // namespace dvecreduce

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
