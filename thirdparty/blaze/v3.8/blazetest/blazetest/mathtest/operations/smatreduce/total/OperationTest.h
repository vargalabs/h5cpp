//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/smatreduce/total/OperationTest.h
//  \brief Header file for the sparse matrix total reduction operation test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_SMATREDUCE_TOTAL_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_SMATREDUCE_TOTAL_OPERATIONTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/DynamicMatrix.h>
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

namespace smatreduce {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the sparse matrix total reduction operation test.
//
// This class template represents one particular test of a total reduction operation on a matrix
// of a particular type. The template argument \a MT represents the type of the matrix operand.
*/
template< typename MT >  // Type of the sparse matrix
class OperationTest
{
 private:
   //**Type definitions****************************************************************************
   using ET = blaze::ElementType_t<MT>;  //!< Element type.

   using OMT  = blaze::OppositeType_t<MT>;    //!< Matrix type with opposite storage order.
   using TMT  = blaze::TransposeType_t<MT>;   //!< Transpose matrix type.
   using TOMT = blaze::TransposeType_t<OMT>;  //!< Transpose matrix type with opposite storage order.

   using RE = blaze::ReduceTrait_t<MT,blaze::Add>;  //!< Result type of the reduction operation.

   using RT = blaze::DynamicMatrix<ET,false>;  //!< Reference type.
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename OP >
   explicit OperationTest( const Creator<MT>& creator, OP op );
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
   MT  mat_;     //!< The sparse matrix operand.
   OMT omat_;    //!< The sparse matrix with opposite storage order.
   RE  res_;     //!< The result of the reduction operation.
   RT  refmat_;  //!< The reference matrix.
   RE  refres_;  //!< The reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OMT  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TMT  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TOMT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( RT   );

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( MT   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OMT  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TMT  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TOMT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( RT   );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET, blaze::ElementType_t<OMT>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET, blaze::ElementType_t<TMT>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET, blaze::ElementType_t<TOMT> );
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
/*!\brief Constructor for the sparse matrix reduction operation test.
//
// \param creator The creator for sparse matrix operand.
// \param op The reduction operation.
// \exception std::runtime_error Operation error detected.
*/
template< typename MT >  // Type of the sparse matrix
template< typename OP >  // Type of the reduction operation
OperationTest<MT>::OperationTest( const Creator<MT>& creator, OP op )
   : mat_( creator() )  // The sparse matrix operand
   , omat_( mat_ )      // The sparse matrix with opposite storage order
   , res_()             // The result of the reduction operation
   , refmat_( mat_ )    // The reference matrix
   , refres_()          // The reference result
   , test_()            // Label of the currently performed test
   , error_()           // Description of the current error type
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
/*!\brief Tests on the initial status of the matrix.
//
// \return void
// \exception std::runtime_error Initialization error detected.
//
// This function runs tests on the initial status of the matrix. In case any initialization
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
void OperationTest<MT>::testInitialStatus()
{
   //=====================================================================================
   // Performing initial tests with the row-major types
   //=====================================================================================

   // Checking the number of rows of the sparse operand
   if( mat_.rows() != refmat_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of row-major sparse operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Detected number of rows = " << mat_.rows() << "\n"
          << "   Expected number of rows = " << refmat_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the sparse operand
   if( mat_.columns() != refmat_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of row-major sparse operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Detected number of columns = " << mat_.columns() << "\n"
          << "   Expected number of columns = " << refmat_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the sparse operand
   if( !isEqual( mat_, refmat_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of row-major sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Current initialization:\n" << mat_ << "\n"
          << "   Expected initialization:\n" << refmat_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing initial tests with the column-major types
   //=====================================================================================

   // Checking the number of rows of the sparse operand
   if( omat_.rows() != refmat_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of column-major sparse operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Detected number of rows = " << omat_.rows() << "\n"
          << "   Expected number of rows = " << refmat_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the sparse operand
   if( omat_.columns() != refmat_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of column-major sparse operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Detected number of columns = " << omat_.columns() << "\n"
          << "   Expected number of columns = " << refmat_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the sparse operand
   if( !isEqual( omat_, refmat_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of column-major sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Current initialization:\n" << omat_ << "\n"
          << "   Expected initialization:\n" << refmat_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the matrix assignment.
//
// \return void
// \exception std::runtime_error Assignment error detected.
//
// This function tests the matrix assignment. In case any error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename MT >  // Type of the sparse matrix
void OperationTest<MT>::testAssignment()
{
   //=====================================================================================
   // Performing an assignment with the row-major types
   //=====================================================================================

   try {
      mat_ = refmat_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the row-major types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( mat_, refmat_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of row-major sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Current initialization:\n" << mat_ << "\n"
          << "   Expected initialization:\n" << refmat_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing an assignment with the column-major types
   //=====================================================================================

   try {
      omat_ = refmat_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the column-major types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( OMT ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( omat_, refmat_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of column-major sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( OMT ).name() << "\n"
          << "   Current initialization:\n" << omat_ << "\n"
          << "   Expected initialization:\n" << refmat_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain sparse matrix reduction operation.
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
template< typename MT >  // Type of the sparse matrix
template< typename OP >  // Type of the reduction operation
void OperationTest<MT>::testBasicOperation( OP op )
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Reduction operation
      //=====================================================================================

      // Reduction operation with the given matrix
      {
         test_  = "Reduction operation with the given matrix";
         error_ = "Failed reduction operation";

         try {
            initResults();
            res_    = reduce( mat_, op );
            refres_ = reduce( refmat_, op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();

         try {
            initResults();
            res_    = reduce( omat_, op );
            refres_ = reduce( refmat_, op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }

      // Reduction operation with evaluated matrix
      {
         test_  = "Reduction operation with evaluated matrices";
         error_ = "Failed reduction operation";

         try {
            initResults();
            res_    = reduce( eval( mat_ ), op );
            refres_ = reduce( eval( refmat_ ), op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();

         try {
            initResults();
            res_    = reduce( eval( omat_ ), op );
            refres_ = reduce( eval( refmat_ ), op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }


      //=====================================================================================
      // Reduction operation with addition assignment
      //=====================================================================================

      // Reduction operation with addition assignment with the given matrix
      {
         test_  = "Reduction operation with addition assignment with the given matrix";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            res_    += reduce( mat_, op );
            refres_ += reduce( refmat_, op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();

         try {
            initResults();
            res_    += reduce( omat_, op );
            refres_ += reduce( refmat_, op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }

      // Reduction operation with addition assignment with evaluated matrix
      {
         test_  = "Reduction operation with addition assignment with evaluated matrix";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            res_    += reduce( eval( mat_ ), op );
            refres_ += reduce( eval( refmat_ ), op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();

         try {
            initResults();
            res_    += reduce( eval( omat_ ), op );
            refres_ += reduce( eval( refmat_ ), op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }


      //=====================================================================================
      // Reduction operation with subtraction assignment
      //=====================================================================================

      // Reduction operation with subtraction assignment with the given matrix
      {
         test_  = "Reduction operation with subtraction assignment with the given matrix";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            res_    -= reduce( mat_, op );
            refres_ -= reduce( refmat_, op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();

         try {
            initResults();
            res_    -= reduce( omat_, op );
            refres_ -= reduce( refmat_, op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }

      // Reduction operation with subtraction assignment with evaluated matrix
      {
         test_  = "Reduction operation with subtraction assignment with evaluated matrix";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            res_    -= reduce( eval( mat_ ), op );
            refres_ -= reduce( eval( refmat_ ), op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();

         try {
            initResults();
            res_    -= reduce( eval( omat_ ), op );
            refres_ -= reduce( eval( refmat_ ), op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }


      //=====================================================================================
      // Reduction operation with multiplication assignment
      //=====================================================================================

      // Reduction operation with multiplication assignment with the given matrix
      {
         test_  = "Reduction operation with multiplication assignment with the given matrix";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            res_    *= reduce( mat_, op );
            refres_ *= reduce( refmat_, op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();

         try {
            initResults();
            res_    *= reduce( omat_, op );
            refres_ *= reduce( refmat_, op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();
      }

      // Reduction operation with multiplication assignment with evaluated matrix
      {
         test_  = "Reduction operation with multiplication assignment with evaluated matrix";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            res_    *= reduce( eval( mat_ ), op );
            refres_ *= reduce( eval( refmat_ ), op );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResult();

         try {
            initResults();
            res_    *= reduce( eval( omat_ ), op );
            refres_ *= reduce( eval( refmat_ ), op );
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
template< typename MT >  // Type of the sparse matrix
void OperationTest<MT>::checkResult()
{
   if( !isEqual( res_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
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
template< typename MT >  // Type of the sparse matrix
void OperationTest<MT>::initResults()
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
template< typename MT >  // Type of the sparse matrix
void OperationTest<MT>::convertException( const std::exception& ex )
{
   std::ostringstream oss;
   oss << " Test : " << test_ << "\n"
       << " Error: " << error_ << "\n"
       << " Details:\n"
       << "   Random seed = " << blaze::getSeed() << "\n"
       << "   Row-major sparse matrix type:\n"
       << "     " << typeid( MT ).name() << "\n"
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
/*!\brief Testing the reduction operation for a specific matrix type.
//
// \param creator The creator for the sparse matrix.
// \return void
*/
template< typename MT >  // Type of the sparse matrix
void runTest( const Creator<MT>& creator )
{
#if BLAZETEST_MATHTEST_TEST_ADDITION
   if( BLAZETEST_MATHTEST_TEST_ADDITION > 1 )
   {
      class Sum : public blaze::Add
      {};

      for( size_t rep=0UL; rep<BLAZETEST_REPETITIONS; ++rep ) {
         OperationTest<MT>( creator, []( const auto& a, const auto& b ){ return a + b; } );
         OperationTest<MT>( creator, blaze::Add() );
         OperationTest<MT>( creator, Sum() );
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
/*!\brief Macro for the definition of a sparse matrix total reduction operation test case.
*/
#define DEFINE_SMATREDUCE_TOTAL_OPERATION_TEST( MT ) \
   extern template class blazetest::mathtest::operations::smatreduce::OperationTest<MT>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a sparse matrix total reduction operation test case.
*/
#define RUN_SMATREDUCE_TOTAL_OPERATION_TEST( C ) \
   blazetest::mathtest::operations::smatreduce::runTest( C )
/*! \endcond */
//*************************************************************************************************

} // namespace smatreduce

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
