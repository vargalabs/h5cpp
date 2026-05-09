//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/dvecrepeat/OperationTest.h
//  \brief Header file for the dense vector repeat operation test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_DVECREPEAT_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_DVECREPEAT_OPERATIONTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>
#include <blaze/math/Aliases.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Functors.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/traits/RepeatTrait.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingScalar.h>
#include <blaze/math/Views.h>
#include <blaze/util/constraints/DerivedFrom.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/Not.h>
#include <blaze/util/Random.h>
#include <blaze/util/typetraits/RemoveCVRef.h>
#include <blazetest/system/MathTest.h>
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/IsEqual.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>


namespace blazetest {

namespace mathtest {

namespace operations {

namespace dvecrepeat {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the dense vector repeat operation test.
//
// This class template represents one particular test of a repeat operation on a vector of a
// particular type. The template argument \a VT represents the type of the vector operand.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
class OperationTest
{
 private:
   //**Type definitions****************************************************************************
   using ET = blaze::ElementType_t<VT>;  //!< Element type.

   using TVT = blaze::TransposeType_t<VT>;  //!< Transpose vector type.

   using DRE   = blaze::RepeatTrait_t<VT,R0>;  //!< Dense result type.
   using DET   = blaze::ElementType_t<DRE>;    //!< Element type of the dense result
   using TDRE  = blaze::TransposeType_t<DRE>;  //!< Transpose dense result type

   using SRE   = blaze::CompressedVector<DET>;  //!< Sparse result type
   using SET   = blaze::ElementType_t<SRE>;     //!< Element type of the sparse result
   using TSRE  = blaze::TransposeType_t<SRE>;   //!< Transpose sparse result type

   using RT  = blaze::CompressedVector<ET>;  //!< Reference type.
   using RRE = blaze::RepeatTrait_t<RT,R0>;  //!< Reference result type.

   using TRT  = blaze::TransposeType_t<RT>;    //!< Transpose reference type.
   using TRRE = blaze::RepeatTrait_t<TRT,R0>;  //!< Transpose reference result type.

   //! Type of the vector repeater expression (runtime argument)
   using VecRepeatExprType1 =
      blaze::RemoveCVRef_t< decltype( blaze::repeat( std::declval<VT>(), R0 ) ) >;

   //! Type of the vector repeater expression (compile time argument)
   using VecRepeatExprType2 =
      blaze::RemoveCVRef_t< decltype( blaze::repeat<R0>( std::declval<VT>() ) ) >;

   //! Type of the transpose vector repeater expression (runtime argument)
   using TVecRepeatExprType1 =
      blaze::RemoveCVRef_t< decltype( blaze::repeat( std::declval<TVT>(), R0 ) ) >;

   //! Type of the transpose vector repeater expression (compile time argument)
   using TVecRepeatExprType2 =
      blaze::RemoveCVRef_t< decltype( blaze::repeat<R0>( std::declval<TVT>() ) ) >;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit OperationTest( const Creator<VT>& creator );
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
                          void testInitialStatus     ();
                          void testAssignment        ();
                          void testEvaluation        ();
                          void testElementAccess     ();
                          void testBasicOperation    ();
                          void testNegatedOperation  ();
   template< typename T > void testScaledOperation   ( T scalar );
                          void testTransOperation    ();
                          void testCTransOperation   ();
                          void testAbsOperation      ();
                          void testConjOperation     ();
                          void testRealOperation     ();
                          void testImagOperation     ();
                          void testEvalOperation     ();
                          void testSerialOperation   ();
                          void testNoAliasOperation  ();
                          void testNoSIMDOperation   ();
                          void testSubvectorOperation( blaze::TrueType  );
                          void testSubvectorOperation( blaze::FalseType );
                          void testElementsOperation ( blaze::TrueType  );
                          void testElementsOperation ( blaze::FalseType );

   template< typename OP > void testCustomOperation( OP op, const std::string& name );
   //@}
   //**********************************************************************************************

   //**Error detection functions*******************************************************************
   /*!\name Error detection functions */
   //@{
   template< typename T > void checkResults();
   template< typename T > void checkTransposeResults();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initResults();
   void initTransposeResults();
   template< typename T > void convertException( const std::exception& ex );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VT    vec_;      //!< The dense vector operand.
   DRE   dres_;     //!< The dense result vector.
   SRE   sres_;     //!< The sparse result vector.
   RT    refvec_;   //!< The reference vector.
   RRE   refres_;   //!< The reference result.
   TVT   tvec_;     //!< The transpose dense vector operand.
   TDRE  tdres_;    //!< The transpose dense result vector.
   TSRE  tsres_;    //!< The transpose sparse result vector.
   TRT   trefvec_;  //!< The transpose reference vector.
   TRRE  trefres_;  //!< The transpose reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TVT  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( RT   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TRT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( DRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( SRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TDRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TSRE );

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( VT    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( TVT   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( RT    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( TRT   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( RRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( TRRE  );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DRE   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TDRE  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SRE   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TSRE  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RRE   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TRRE  );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET , blaze::ElementType_t<TVT>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<DRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<TDRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<SRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<SRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<TSRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<DRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT , blaze::TransposeType_t<TVT> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RT , blaze::TransposeType_t<TRT> );

   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG     ( VecRepeatExprType1, blaze::ResultType_t<VecRepeatExprType1>    );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_DIFFERENT_TRANSPOSE_FLAG( VecRepeatExprType1, blaze::TransposeType_t<VecRepeatExprType1> );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG     ( VecRepeatExprType2, blaze::ResultType_t<VecRepeatExprType2>    );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_DIFFERENT_TRANSPOSE_FLAG( VecRepeatExprType2, blaze::TransposeType_t<VecRepeatExprType2> );

   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG     ( TVecRepeatExprType1, blaze::ResultType_t<TVecRepeatExprType1>    );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_DIFFERENT_TRANSPOSE_FLAG( TVecRepeatExprType1, blaze::TransposeType_t<TVecRepeatExprType1> );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG     ( TVecRepeatExprType2, blaze::ResultType_t<TVecRepeatExprType2>    );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_DIFFERENT_TRANSPOSE_FLAG( TVecRepeatExprType2, blaze::TransposeType_t<TVecRepeatExprType2> );

   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( VecRepeatExprType1 , blaze::BaseType_t<VecRepeatExprType1>  );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( VecRepeatExprType2 , blaze::BaseType_t<VecRepeatExprType2>  );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TVecRepeatExprType1, blaze::BaseType_t<TVecRepeatExprType1> );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TVecRepeatExprType2, blaze::BaseType_t<TVecRepeatExprType2> );
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
/*!\brief Constructor for the dense vector repeat operation test.
//
// \param creator The creator for dense vector operand.
// \exception std::runtime_error Operation error detected.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
OperationTest<VT,R0>::OperationTest( const Creator<VT>& creator )
   : vec_( creator() )     // The dense vector operand
   , dres_()               // The dense result vector
   , sres_()               // The sparse result vector
   , refvec_( vec_ )       // The reference vector
   , refres_()             // The reference result
   , tvec_( trans(vec_) )  // The transpose dense vector operand
   , tdres_()              // The transpose dense result vector
   , tsres_()              // The transpose sparse result vector
   , trefvec_( tvec_ )     // The transpose reference vector
   , trefres_()            // The transpose reference result
   , test_()               // Label of the currently performed test
   , error_()              // Description of the current error type
{
   using namespace blaze;

   using Scalar = UnderlyingScalar_t<DET>;

   testInitialStatus();
   testAssignment();
   testEvaluation();
   testElementAccess();
   testBasicOperation();
   testNegatedOperation();
   testScaledOperation( 2 );
   testScaledOperation( 2UL );
   testScaledOperation( 2.0F );
   testScaledOperation( 2.0 );
   testScaledOperation( Scalar( 2 ) );
   testTransOperation();
   testCTransOperation();
   testAbsOperation();
   testConjOperation();
   testRealOperation();
   testImagOperation();
   testEvalOperation();
   testSerialOperation();
   testNoAliasOperation();
   testNoSIMDOperation();
   testSubvectorOperation( Not_t< IsUniform<DRE> >() );
   testElementsOperation( Not_t< IsUniform<DRE> >() );
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Tests on the initial status of the vector.
//
// \return void
// \exception std::runtime_error Initialization error detected.
//
// This function runs tests on the initial status of the vector. In case any initialization
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testInitialStatus()
{
   //=====================================================================================
   // Performing initial tests with the given vector
   //=====================================================================================

   // Checking the size of the vector operand
   if( vec_.size() != refvec_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of dense vector operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Detected size = " << vec_.size() << "\n"
          << "   Expected size = " << refvec_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the vector operand
   if( !isEqual( vec_, refvec_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of dense vector operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Current initialization:\n" << vec_ << "\n"
          << "   Expected initialization:\n" << refvec_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing initial tests with the transpose type
   //=====================================================================================

   // Checking the size of the vector operand
   if( tvec_.size() != trefvec_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of transpose dense vector operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT ).name() << "\n"
          << "   Detected size = " << tvec_.size() << "\n"
          << "   Expected size = " << trefvec_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the vector operand
   if( !isEqual( tvec_, trefvec_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of transpose dense vector operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT ).name() << "\n"
          << "   Current initialization:\n" << tvec_ << "\n"
          << "   Expected initialization:\n" << trefvec_ << "\n";
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
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testAssignment()
{
   //=====================================================================================
   // Performing an assignment with the given vector
   //=====================================================================================

   try {
      vec_ = refvec_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the given vectors\n"
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
      oss << " Test: Checking the assignment result of dense vector operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Current initialization:\n" << vec_ << "\n"
          << "   Expected initialization:\n" << refvec_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing an assignment with the transpose type
   //=====================================================================================

   try {
      tvec_ = trefvec_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the transpose types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( vec_, refvec_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of transpose dense vector operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT ).name() << "\n"
          << "   Current initialization:\n" << tvec_ << "\n"
          << "   Expected initialization:\n" << trefvec_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the explicit evaluation.
//
// \return void
// \exception std::runtime_error Evaluation error detected.
//
// This function tests the explicit evaluation. In case any error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testEvaluation()
{
   using blaze::repeat;


   //=====================================================================================
   // Testing the evaluation with a column vector
   //=====================================================================================

   {
      const auto res   ( evaluate( repeat( vec_, R0 ) ) );
      const auto refres( evaluate( repeat( refvec_, R0 ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given vector (runtime)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( vec_ ).name() << "\n"
             << "   Deduced result type:\n"
             << "     " << typeid( res ).name() << "\n"
             << "   Deduced reference result type:\n"
             << "     " << typeid( refres ).name() << "\n"
             << "   Result:\n" << res << "\n"
             << "   Expected result:\n" << refres << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      const auto res   ( evaluate( repeat<R0>( vec_ ) ) );
      const auto refres( evaluate( repeat<R0>( refvec_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given vector (compile time)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( vec_ ).name() << "\n"
             << "   Deduced result type:\n"
             << "     " << typeid( res ).name() << "\n"
             << "   Deduced reference result type:\n"
             << "     " << typeid( refres ).name() << "\n"
             << "   Result:\n" << res << "\n"
             << "   Expected result:\n" << refres << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      const auto res   ( evaluate( repeat( eval( vec_ ), R0 ) ) );
      const auto refres( evaluate( repeat( eval( refvec_ ), R0 ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated vector (runtime)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( vec_ ).name() << "\n"
             << "   Deduced result type:\n"
             << "     " << typeid( res ).name() << "\n"
             << "   Deduced reference result type:\n"
             << "     " << typeid( refres ).name() << "\n"
             << "   Result:\n" << res << "\n"
             << "   Expected result:\n" << refres << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      const auto res   ( evaluate( repeat<R0>( eval( vec_ ) ) ) );
      const auto refres( evaluate( repeat<R0>( eval( refvec_ ) ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated vector (compile time)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( vec_ ).name() << "\n"
             << "   Deduced result type:\n"
             << "     " << typeid( res ).name() << "\n"
             << "   Deduced reference result type:\n"
             << "     " << typeid( refres ).name() << "\n"
             << "   Result:\n" << res << "\n"
             << "   Expected result:\n" << refres << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Testing the evaluation with a row vector
   //=====================================================================================

   {
      const auto res   ( evaluate( repeat( tvec_, R0 ) ) );
      const auto refres( evaluate( repeat( trefvec_, R0 ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given vector (runtime)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense row vector type:\n"
             << "     " << typeid( vec_ ).name() << "\n"
             << "   Deduced result type:\n"
             << "     " << typeid( res ).name() << "\n"
             << "   Deduced reference result type:\n"
             << "     " << typeid( refres ).name() << "\n"
             << "   Result:\n" << res << "\n"
             << "   Expected result:\n" << refres << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      const auto res   ( evaluate( repeat<R0>( tvec_ ) ) );
      const auto refres( evaluate( repeat<R0>( trefvec_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given vector (compile time)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense row vector type:\n"
             << "     " << typeid( vec_ ).name() << "\n"
             << "   Deduced result type:\n"
             << "     " << typeid( res ).name() << "\n"
             << "   Deduced reference result type:\n"
             << "     " << typeid( refres ).name() << "\n"
             << "   Result:\n" << res << "\n"
             << "   Expected result:\n" << refres << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      const auto res   ( evaluate( repeat( eval( tvec_ ), R0 ) ) );
      const auto refres( evaluate( repeat( eval( trefvec_ ), R0 ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated vector (runtime)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense row vector type:\n"
             << "     " << typeid( vec_ ).name() << "\n"
             << "   Deduced result type:\n"
             << "     " << typeid( res ).name() << "\n"
             << "   Deduced reference result type:\n"
             << "     " << typeid( refres ).name() << "\n"
             << "   Result:\n" << res << "\n"
             << "   Expected result:\n" << refres << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      const auto res   ( evaluate( repeat<R0>( eval( tvec_ ) ) ) );
      const auto refres( evaluate( repeat<R0>( eval( trefvec_ ) ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated vector (compile time)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense row vector type:\n"
             << "     " << typeid( vec_ ).name() << "\n"
             << "   Deduced result type:\n"
             << "     " << typeid( res ).name() << "\n"
             << "   Deduced reference result type:\n"
             << "     " << typeid( refres ).name() << "\n"
             << "   Result:\n" << res << "\n"
             << "   Expected result:\n" << refres << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the vector element access.
//
// \return void
// \exception std::runtime_error Element access error detected.
//
// This function tests the element access via the subscript operator. In case any
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testElementAccess()
{
   using blaze::equal;
   using blaze::repeat;


   //=====================================================================================
   // Testing the element access with a column vector
   //=====================================================================================

   if( vec_.size() > 0UL && R0 > 0UL )
   {
      const size_t n( vec_.size()*R0 - 1UL );

      if( !equal( repeat( vec_, R0 )[n], repeat( refvec_, R0 )[n] ) ||
          !equal( repeat( vec_, R0 ).at(n), repeat( refvec_, R0 ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of repeater expression (runtime)\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( repeat<R0>( vec_ )[n], repeat<R0>( refvec_ )[n] ) ||
          !equal( repeat<R0>( vec_ ).at(n), repeat<R0>( refvec_ ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of repeater expression (compile time)\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( repeat( eval( vec_ ), R0 )[n], repeat( refvec_, R0 )[n] ) ||
          !equal( repeat( eval( vec_ ), R0 ).at(n), repeat( refvec_, R0 ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of evaluated repeater expression (runtime)\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( repeat<R0>( eval( vec_ ) )[n], repeat<R0>( refvec_ )[n] ) ||
          !equal( repeat<R0>( eval( vec_ ) ).at(n), repeat<R0>( refvec_ ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of evaluated repeater expression (compile time)\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      repeat( vec_, R0 ).at( vec_.size()*R0 );

      std::ostringstream oss;
      oss << " Test : Checked element access of repeater expression (runtime)\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense column vector type:\n"
          << "     " << typeid( VT ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      repeat<R0>( vec_ ).at( vec_.size()*R0 );

      std::ostringstream oss;
      oss << " Test : Checked element access of repeater expression (compile time)\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense column vector type:\n"
          << "     " << typeid( VT ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}


   //=====================================================================================
   // Testing the element access with a row vector
   //=====================================================================================

   if( tvec_.size() > 0UL && R0 > 0UL )
   {
      const size_t n( tvec_.size()*R0 - 1UL );

      if( !equal( repeat( tvec_, R0 )[n], repeat( trefvec_, R0 )[n] ) ||
          !equal( repeat( tvec_, R0 ).at(n), repeat( trefvec_, R0 ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of repeater expression (runtime)\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( TVT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( repeat<R0>( tvec_ )[n], repeat<R0>( trefvec_ )[n] ) ||
          !equal( repeat<R0>( tvec_ ).at(n), repeat<R0>( trefvec_ ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of repeater expression (compile time)\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( TVT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( repeat( eval( tvec_ ), R0 )[n], repeat( trefvec_, R0 )[n] ) ||
          !equal( repeat( eval( tvec_ ), R0 ).at(n), repeat( trefvec_, R0 ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of evaluated repeater expression (runtime)\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( TVT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( repeat<R0>( eval( tvec_ ) )[n], repeat<R0>( trefvec_ )[n] ) ||
          !equal( repeat<R0>( eval( tvec_ ) ).at(n), repeat<R0>( trefvec_ ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of evaluated repeater expression (compile time)\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( TVT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      repeat( tvec_, R0 ).at( tvec_.size()*R0 );

      std::ostringstream oss;
      oss << " Test : Checked element access of repeater expression\n"
          << " Error: Out-of-bound access succeeded (runtime)\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense column vector type:\n"
          << "     " << typeid( TVT ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      repeat<R0>( tvec_ ).at( tvec_.size()*R0 );

      std::ostringstream oss;
      oss << " Test : Checked element access of repeater expression\n"
          << " Error: Out-of-bound access succeeded (compile time)\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense column vector type:\n"
          << "     " << typeid( TVT ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the vector repeat operation with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case
// any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      using blaze::repeat;


      //=====================================================================================
      // Repeat operation
      //=====================================================================================

      // Repeat operation with the given vector (runtime)
      {
         test_  = "Repeat operation with the given vector (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat( vec_, R0 );
            sres_   = repeat( vec_, R0 );
            refres_ = repeat( refvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = repeat( tvec_, R0 );
            tsres_   = repeat( tvec_, R0 );
            trefres_ = repeat( trefvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Repeat operation with the given vector (compile time)
      {
         test_  = "Repeat operation with the given vector (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat<R0>( vec_ );
            sres_   = repeat<R0>( vec_ );
            refres_ = repeat<R0>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = repeat<R0>( tvec_ );
            tsres_   = repeat<R0>( tvec_ );
            trefres_ = repeat<R0>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Repeat operation with evaluated vector (runtime)
      {
         test_  = "Repeat operation with evaluated vector (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat( eval( vec_ ), R0 );
            sres_   = repeat( eval( vec_ ), R0 );
            refres_ = repeat( eval( refvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = repeat( eval( tvec_ ), R0 );
            tsres_   = repeat( eval( tvec_ ), R0 );
            trefres_ = repeat( eval( trefvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Repeat operation with evaluated vector (compile time)
      {
         test_  = "Repeat operation with evaluated vector (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat<R0>( eval( vec_ ) );
            sres_   = repeat<R0>( eval( vec_ ) );
            refres_ = repeat<R0>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = repeat<R0>( eval( tvec_ ) );
            tsres_   = repeat<R0>( eval( tvec_ ) );
            trefres_ = repeat<R0>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Repeat with addition assignment
      //=====================================================================================

      // Repeat with addition assignment with the given vector (runtime)
      {
         test_  = "Repeat with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat( vec_, R0 );
            sres_   += repeat( vec_, R0 );
            refres_ += repeat( refvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += repeat( tvec_, R0 );
            tsres_   += repeat( tvec_, R0 );
            trefres_ += repeat( trefvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Repeat with addition assignment with the given vector (compile time)
      {
         test_  = "Repeat with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat<R0>( vec_ );
            sres_   += repeat<R0>( vec_ );
            refres_ += repeat<R0>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += repeat<R0>( tvec_ );
            tsres_   += repeat<R0>( tvec_ );
            trefres_ += repeat<R0>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Repeat with addition assignment with evaluated vector (runtime)
      {
         test_  = "Repeat with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat( eval( vec_ ), R0 );
            sres_   += repeat( eval( vec_ ), R0 );
            refres_ += repeat( eval( refvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += repeat( eval( tvec_ ), R0 );
            tsres_   += repeat( eval( tvec_ ), R0 );
            trefres_ += repeat( eval( trefvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Repeat with addition assignment with evaluated vector (compile time)
      {
         test_  = "Repeat with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat<R0>( eval( vec_ ) );
            sres_   += repeat<R0>( eval( vec_ ) );
            refres_ += repeat<R0>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += repeat<R0>( eval( tvec_ ) );
            tsres_   += repeat<R0>( eval( tvec_ ) );
            trefres_ += repeat<R0>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Repeat with subtraction assignment
      //=====================================================================================

      // Repeat with subtraction assignment with the given vector (runtime)
      {
         test_  = "Repeat with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat( vec_, R0 );
            sres_   -= repeat( vec_, R0 );
            refres_ -= repeat( refvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= repeat( tvec_, R0 );
            tsres_   -= repeat( tvec_, R0 );
            trefres_ -= repeat( trefvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Repeat with subtraction assignment with the given vector (compile time)
      {
         test_  = "Repeat with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat<R0>( vec_ );
            sres_   -= repeat<R0>( vec_ );
            refres_ -= repeat<R0>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= repeat<R0>( tvec_ );
            tsres_   -= repeat<R0>( tvec_ );
            trefres_ -= repeat<R0>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Repeat with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Repeat with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat( eval( vec_ ), R0 );
            sres_   -= repeat( eval( vec_ ), R0 );
            refres_ -= repeat( eval( refvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= repeat( eval( tvec_ ), R0 );
            tsres_   -= repeat( eval( tvec_ ), R0 );
            trefres_ -= repeat( eval( trefvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Repeat with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Repeat with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat<R0>( eval( vec_ ) );
            sres_   -= repeat<R0>( eval( vec_ ) );
            refres_ -= repeat<R0>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= repeat<R0>( eval( tvec_ ) );
            tsres_   -= repeat<R0>( eval( tvec_ ) );
            trefres_ -= repeat<R0>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Repeat with multiplication assignment
      //=====================================================================================

      // Repeat with multiplication assignment with the given vector (runtime)
      {
         test_  = "Repeat with multiplication assignment with the given vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= repeat( vec_, R0 );
            sres_   *= repeat( vec_, R0 );
            refres_ *= repeat( refvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= repeat( tvec_, R0 );
            tsres_   *= repeat( tvec_, R0 );
            trefres_ *= repeat( trefvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Repeat with multiplication assignment with the given vector (compile time)
      {
         test_  = "Repeat with multiplication assignment with the given vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= repeat<R0>( vec_ );
            sres_   *= repeat<R0>( vec_ );
            refres_ *= repeat<R0>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= repeat<R0>( tvec_ );
            tsres_   *= repeat<R0>( tvec_ );
            trefres_ *= repeat<R0>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Repeat with multiplication assignment with evaluated vector (runtime)
      {
         test_  = "Repeat with multiplication assignment with evaluated vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= repeat( eval( vec_ ), R0 );
            sres_   *= repeat( eval( vec_ ), R0 );
            refres_ *= repeat( eval( refvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= repeat( eval( tvec_ ), R0 );
            tsres_   *= repeat( eval( tvec_ ), R0 );
            trefres_ *= repeat( eval( trefvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Repeat with multiplication assignment with evaluated vector (compile time)
      {
         test_  = "Repeat with multiplication assignment with evaluated vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= repeat<R0>( eval( vec_ ) );
            sres_   *= repeat<R0>( eval( vec_ ) );
            refres_ *= repeat<R0>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= repeat<R0>( eval( tvec_ ) );
            tsres_   *= repeat<R0>( eval( tvec_ ) );
            trefres_ *= repeat<R0>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the negated dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the negated vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment, and division assignment.
// In case any error resulting from the repeat operation or the subsequent assignment is
// detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      using blaze::repeat;


      //=====================================================================================
      // Negated repeat operation
      //=====================================================================================

      // Negated repeat operation with the given vector (runtime)
      {
         test_  = "Negated repeat operation with the given vector (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = -repeat( vec_, R0 );
            sres_   = -repeat( vec_, R0 );
            refres_ = -repeat( refvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = -repeat( tvec_, R0 );
            tsres_   = -repeat( tvec_, R0 );
            trefres_ = -repeat( trefvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated repeat operation with the given vector (compile time)
      {
         test_  = "Negated repeat operation with the given vector (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = -repeat<R0>( vec_ );
            sres_   = -repeat<R0>( vec_ );
            refres_ = -repeat<R0>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = -repeat<R0>( tvec_ );
            tsres_   = -repeat<R0>( tvec_ );
            trefres_ = -repeat<R0>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated repeat operation with evaluated vector (runtime)
      {
         test_  = "Negated repeat operation with evaluated vector (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = -repeat( eval( vec_ ), R0 );
            sres_   = -repeat( eval( vec_ ), R0 );
            refres_ = -repeat( eval( refvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = -repeat( eval( tvec_ ), R0 );
            tsres_   = -repeat( eval( tvec_ ), R0 );
            trefres_ = -repeat( eval( trefvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated repeat operation with evaluated vector (compile time)
      {
         test_  = "Negated repeat operation with evaluated vector (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = -repeat<R0>( eval( vec_ ) );
            sres_   = -repeat<R0>( eval( vec_ ) );
            refres_ = -repeat<R0>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = -repeat<R0>( eval( tvec_ ) );
            tsres_   = -repeat<R0>( eval( tvec_ ) );
            trefres_ = -repeat<R0>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Negated repeat with addition assignment
      //=====================================================================================

      // Negated repeat with addition assignment with the given vector (runtime)
      {
         test_  = "Negated repeat with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += -repeat( vec_, R0 );
            sres_   += -repeat( vec_, R0 );
            refres_ += -repeat( refvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += -repeat( tvec_, R0 );
            tsres_   += -repeat( tvec_, R0 );
            trefres_ += -repeat( trefvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated repeat with addition assignment with the given vector (compile time)
      {
         test_  = "Negated repeat with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += -repeat<R0>( vec_ );
            sres_   += -repeat<R0>( vec_ );
            refres_ += -repeat<R0>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += -repeat<R0>( tvec_ );
            tsres_   += -repeat<R0>( tvec_ );
            trefres_ += -repeat<R0>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated repeat with addition assignment with evaluated vector (runtime)
      {
         test_  = "Negated repeat with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += -repeat( eval( vec_ ), R0 );
            sres_   += -repeat( eval( vec_ ), R0 );
            refres_ += -repeat( eval( refvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += -repeat( eval( tvec_ ), R0 );
            tsres_   += -repeat( eval( tvec_ ), R0 );
            trefres_ += -repeat( eval( trefvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated repeat with addition assignment with evaluated vector (compile time)
      {
         test_  = "Negated repeat with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += -repeat<R0>( eval( vec_ ) );
            sres_   += -repeat<R0>( eval( vec_ ) );
            refres_ += -repeat<R0>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += -repeat<R0>( eval( tvec_ ) );
            tsres_   += -repeat<R0>( eval( tvec_ ) );
            trefres_ += -repeat<R0>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Negated repeat with subtraction assignment
      //=====================================================================================

      // Negated repeat with subtraction assignment with the given vector (runtime)
      {
         test_  = "Negated repeat with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= -repeat( vec_, R0 );
            sres_   -= -repeat( vec_, R0 );
            refres_ -= -repeat( refvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= -repeat( tvec_, R0 );
            tsres_   -= -repeat( tvec_, R0 );
            trefres_ -= -repeat( trefvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated repeat with subtraction assignment with the given vector (compile time)
      {
         test_  = "Negated repeat with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= -repeat<R0>( vec_ );
            sres_   -= -repeat<R0>( vec_ );
            refres_ -= -repeat<R0>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= -repeat<R0>( tvec_ );
            tsres_   -= -repeat<R0>( tvec_ );
            trefres_ -= -repeat<R0>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated repeat with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Negated repeat with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= -repeat( eval( vec_ ), R0 );
            sres_   -= -repeat( eval( vec_ ), R0 );
            refres_ -= -repeat( eval( refvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= -repeat( eval( tvec_ ), R0 );
            tsres_   -= -repeat( eval( tvec_ ), R0 );
            trefres_ -= -repeat( eval( trefvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated repeat with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Negated repeat with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= -repeat<R0>( eval( vec_ ) );
            sres_   -= -repeat<R0>( eval( vec_ ) );
            refres_ -= -repeat<R0>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= -repeat<R0>( eval( tvec_ ) );
            tsres_   -= -repeat<R0>( eval( tvec_ ) );
            trefres_ -= -repeat<R0>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Negated repeat with multiplication assignment
      //=====================================================================================

      // Negated repeat with multiplication assignment with the given vector (runtime)
      {
         test_  = "Negated repeat with multiplication assignment with the given vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= -repeat( vec_, R0 );
            sres_   *= -repeat( vec_, R0 );
            refres_ *= -repeat( refvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= -repeat( tvec_, R0 );
            tsres_   *= -repeat( tvec_, R0 );
            trefres_ *= -repeat( trefvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated repeat with multiplication assignment with the given vector (compile time)
      {
         test_  = "Negated repeat with multiplication assignment with the given vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= -repeat<R0>( vec_ );
            sres_   *= -repeat<R0>( vec_ );
            refres_ *= -repeat<R0>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= -repeat<R0>( tvec_ );
            tsres_   *= -repeat<R0>( tvec_ );
            trefres_ *= -repeat<R0>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated repeat with multiplication assignment with evaluated vector (runtime)
      {
         test_  = "Negated repeat with multiplication assignment with evaluated vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= -repeat( eval( vec_ ), R0 );
            sres_   *= -repeat( eval( vec_ ), R0 );
            refres_ *= -repeat( eval( refvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= -repeat( eval( tvec_ ), R0 );
            tsres_   *= -repeat( eval( tvec_ ), R0 );
            trefres_ *= -repeat( eval( trefvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated repeat with multiplication assignment with evaluated vector (compile time)
      {
         test_  = "Negated repeat with multiplication assignment with evaluated vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= -repeat<R0>( eval( vec_ ) );
            sres_   *= -repeat<R0>( eval( vec_ ) );
            refres_ *= -repeat<R0>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= -repeat<R0>( eval( tvec_ ) );
            tsres_   *= -repeat<R0>( eval( tvec_ ) );
            trefres_ *= -repeat<R0>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the scaled dense vector repeat operation.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the scaled vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
template< typename T >    // Type of the scalar
void OperationTest<VT,R0>::testScaledOperation( T scalar )
{
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE( T );

   if( scalar == T(0) )
      throw std::invalid_argument( "Invalid scalar parameter" );


#if BLAZETEST_MATHTEST_TEST_SCALED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SCALED_OPERATION > 1 )
   {
      using blaze::repeat;


      //=====================================================================================
      // Scaled repeat operation (s*OP)
      //=====================================================================================

      // Scaled repeat operation with the given vector (s*OP, runtime)
      {
         test_  = "Scaled repeat operation with the given vector (s*OP, runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = scalar * repeat( vec_, R0 );
            sres_   = scalar * repeat( vec_, R0 );
            refres_ = scalar * repeat( refvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = scalar * repeat( tvec_, R0 );
            tsres_   = scalar * repeat( tvec_, R0 );
            trefres_ = scalar * repeat( trefvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat operation with the given vector (s*OP, compile time)
      {
         test_  = "Scaled repeat operation with the given vector (s*OP, compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = scalar * repeat<R0>( vec_ );
            sres_   = scalar * repeat<R0>( vec_ );
            refres_ = scalar * repeat<R0>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = scalar * repeat<R0>( tvec_ );
            tsres_   = scalar * repeat<R0>( tvec_ );
            trefres_ = scalar * repeat<R0>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat operation with evaluated vector (s*OP, runtime)
      {
         test_  = "Scaled repeat operation with evaluated vector (s*OP, runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = scalar * repeat( eval( vec_ ), R0 );
            sres_   = scalar * repeat( eval( vec_ ), R0 );
            refres_ = scalar * repeat( eval( refvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = scalar * repeat( eval( tvec_ ), R0 );
            tsres_   = scalar * repeat( eval( tvec_ ), R0 );
            trefres_ = scalar * repeat( eval( trefvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat operation with evaluated vector (s*OP, compile time)
      {
         test_  = "Scaled repeat operation with evaluated vector (s*OP, compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = scalar * repeat<R0>( eval( vec_ ) );
            sres_   = scalar * repeat<R0>( eval( vec_ ) );
            refres_ = scalar * repeat<R0>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = scalar * repeat<R0>( eval( tvec_ ) );
            tsres_   = scalar * repeat<R0>( eval( tvec_ ) );
            trefres_ = scalar * repeat<R0>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled repeat operation (OP*s)
      //=====================================================================================

      // Scaled repeat operation with the given vector (OP*s, runtime)
      {
         test_  = "Scaled repeat operation with the given vector (OP*s, runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat( vec_, R0 ) * scalar;
            sres_   = repeat( vec_, R0 ) * scalar;
            refres_ = repeat( refvec_, R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = repeat( tvec_, R0 ) * scalar;
            tsres_   = repeat( tvec_, R0 ) * scalar;
            trefres_ = repeat( trefvec_, R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat operation with the given vector (OP*s, compile time)
      {
         test_  = "Scaled repeat operation with the given vector (OP*s, compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat<R0>( vec_ ) * scalar;
            sres_   = repeat<R0>( vec_ ) * scalar;
            refres_ = repeat<R0>( refvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = repeat<R0>( tvec_ ) * scalar;
            tsres_   = repeat<R0>( tvec_ ) * scalar;
            trefres_ = repeat<R0>( trefvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat operation with evaluated vector (OP*s, runtime)
      {
         test_  = "Scaled repeat operation with evaluated vector (OP*s, runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat( eval( vec_ ), R0 ) * scalar;
            sres_   = repeat( eval( vec_ ), R0 ) * scalar;
            refres_ = repeat( eval( refvec_ ), R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = repeat( eval( tvec_ ), R0 ) * scalar;
            tsres_   = repeat( eval( tvec_ ), R0 ) * scalar;
            trefres_ = repeat( eval( trefvec_ ), R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat operation with evaluated vector (OP*s, compile time)
      {
         test_  = "Scaled repeat operation with evaluated vector (OP*s, compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat<R0>( eval( vec_ ) ) * scalar;
            sres_   = repeat<R0>( eval( vec_ ) ) * scalar;
            refres_ = repeat<R0>( eval( refvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = repeat<R0>( eval( tvec_ ) ) * scalar;
            tsres_   = repeat<R0>( eval( tvec_ ) ) * scalar;
            trefres_ = repeat<R0>( eval( trefvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled repeat operation (OP/s)
      //=====================================================================================

      // Scaled repeat operation with the given vector (OP/s, runtime)
      {
         test_  = "Scaled repeat operation with the given vector (OP/s, runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat( vec_, R0 ) / scalar;
            sres_   = repeat( vec_, R0 ) / scalar;
            refres_ = repeat( refvec_, R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = repeat( tvec_, R0 ) / scalar;
            tsres_   = repeat( tvec_, R0 ) / scalar;
            trefres_ = repeat( trefvec_, R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat operation with the given vector (OP/s, compile time)
      {
         test_  = "Scaled repeat operation with the given vector (OP/s, compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat<R0>( vec_ ) / scalar;
            sres_   = repeat<R0>( vec_ ) / scalar;
            refres_ = repeat<R0>( refvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = repeat<R0>( tvec_ ) / scalar;
            tsres_   = repeat<R0>( tvec_ ) / scalar;
            trefres_ = repeat<R0>( trefvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat operation with evaluated vector (OP/s, runtime)
      {
         test_  = "Scaled repeat operation with evaluated vector (OP/s, runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat( eval( vec_ ), R0 ) / scalar;
            sres_   = repeat( eval( vec_ ), R0 ) / scalar;
            refres_ = repeat( eval( refvec_ ), R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = repeat( eval( tvec_ ), R0 ) / scalar;
            tsres_   = repeat( eval( tvec_ ), R0 ) / scalar;
            trefres_ = repeat( eval( trefvec_ ), R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat operation with evaluated vector (OP/s, compile time)
      {
         test_  = "Scaled repeat operation with evaluated vector (OP/s, compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat<R0>( eval( vec_ ) ) / scalar;
            sres_   = repeat<R0>( eval( vec_ ) ) / scalar;
            refres_ = repeat<R0>( eval( refvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = repeat<R0>( eval( tvec_ ) ) / scalar;
            tsres_   = repeat<R0>( eval( tvec_ ) ) / scalar;
            trefres_ = repeat<R0>( eval( trefvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled repeat with addition assignment (s*OP)
      //=====================================================================================

      // Scaled repeat with addition assignment with the given vector (s*OP, runtime)
      {
         test_  = "Scaled repeat with addition assignment with the given vector (s*OP, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += scalar * repeat( vec_, R0 );
            sres_   += scalar * repeat( vec_, R0 );
            refres_ += scalar * repeat( refvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += scalar * repeat( tvec_, R0 );
            tsres_   += scalar * repeat( tvec_, R0 );
            trefres_ += scalar * repeat( trefvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with addition assignment with the given vector (s*OP, compile time)
      {
         test_  = "Scaled repeat with addition assignment with the given vector (s*OP, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += scalar * repeat<R0>( vec_ );
            sres_   += scalar * repeat<R0>( vec_ );
            refres_ += scalar * repeat<R0>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += scalar * repeat<R0>( tvec_ );
            tsres_   += scalar * repeat<R0>( tvec_ );
            trefres_ += scalar * repeat<R0>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with addition assignment with evaluated vector (s*OP, runtime)
      {
         test_  = "Scaled repeat with addition assignment with evaluated vector (s*OP, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += scalar * repeat( eval( vec_ ), R0 );
            sres_   += scalar * repeat( eval( vec_ ), R0 );
            refres_ += scalar * repeat( eval( refvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += scalar * repeat( eval( tvec_ ), R0 );
            tsres_   += scalar * repeat( eval( tvec_ ), R0 );
            trefres_ += scalar * repeat( eval( trefvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with addition assignment with evaluated vector (s*OP, compile time)
      {
         test_  = "Scaled repeat with addition assignment with evaluated vector (s*OP, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += scalar * repeat<R0>( eval( vec_ ) );
            sres_   += scalar * repeat<R0>( eval( vec_ ) );
            refres_ += scalar * repeat<R0>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += scalar * repeat<R0>( eval( tvec_ ) );
            tsres_   += scalar * repeat<R0>( eval( tvec_ ) );
            trefres_ += scalar * repeat<R0>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled repeat with addition assignment (OP*s)
      //=====================================================================================

      // Scaled repeat with addition assignment with the given vector (OP*s, runtime)
      {
         test_  = "Scaled repeat with addition assignment with the given vector (OP*s, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat( vec_, R0 ) * scalar;
            sres_   += repeat( vec_, R0 ) * scalar;
            refres_ += repeat( refvec_, R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += repeat( tvec_, R0 ) * scalar;
            tsres_   += repeat( tvec_, R0 ) * scalar;
            trefres_ += repeat( trefvec_, R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with addition assignment with the given vector (OP*s, compile time)
      {
         test_  = "Scaled repeat with addition assignment with the given vector (OP*s, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat<R0>( vec_ ) * scalar;
            sres_   += repeat<R0>( vec_ ) * scalar;
            refres_ += repeat<R0>( refvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += repeat<R0>( tvec_ ) * scalar;
            tsres_   += repeat<R0>( tvec_ ) * scalar;
            trefres_ += repeat<R0>( trefvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with addition assignment with evaluated vector (OP*s, runtime)
      {
         test_  = "Scaled repeat with addition assignment with evaluated vector (OP*s, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat( eval( vec_ ), R0 ) * scalar;
            sres_   += repeat( eval( vec_ ), R0 ) * scalar;
            refres_ += repeat( eval( refvec_ ), R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += repeat( eval( tvec_ ), R0 ) * scalar;
            tsres_   += repeat( eval( tvec_ ), R0 ) * scalar;
            trefres_ += repeat( eval( trefvec_ ), R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with addition assignment with evaluated vector (OP*s, compile time)
      {
         test_  = "Scaled repeat with addition assignment with evaluated vector (OP*s, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat<R0>( eval( vec_ ) ) * scalar;
            sres_   += repeat<R0>( eval( vec_ ) ) * scalar;
            refres_ += repeat<R0>( eval( refvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += repeat<R0>( eval( tvec_ ) ) * scalar;
            tsres_   += repeat<R0>( eval( tvec_ ) ) * scalar;
            trefres_ += repeat<R0>( eval( trefvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled repeat with addition assignment (OP/s)
      //=====================================================================================

      // Scaled repeat with addition assignment with the given vector (OP/s, runtime)
      {
         test_  = "Scaled repeat with addition assignment with the given vector (OP/s, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat( vec_, R0 ) / scalar;
            sres_   += repeat( vec_, R0 ) / scalar;
            refres_ += repeat( refvec_, R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += repeat( tvec_, R0 ) / scalar;
            tsres_   += repeat( tvec_, R0 ) / scalar;
            trefres_ += repeat( trefvec_, R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with addition assignment with the given vector (OP/s, compile time)
      {
         test_  = "Scaled repeat with addition assignment with the given vector (OP/s, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat<R0>( vec_ ) / scalar;
            sres_   += repeat<R0>( vec_ ) / scalar;
            refres_ += repeat<R0>( refvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += repeat<R0>( tvec_ ) / scalar;
            tsres_   += repeat<R0>( tvec_ ) / scalar;
            trefres_ += repeat<R0>( trefvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with addition assignment with evaluated vector (OP/s, runtime)
      {
         test_  = "Scaled repeat with addition assignment with evaluated vector (OP/s, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat( eval( vec_ ), R0 ) / scalar;
            sres_   += repeat( eval( vec_ ), R0 ) / scalar;
            refres_ += repeat( eval( refvec_ ), R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += repeat( eval( tvec_ ), R0 ) / scalar;
            tsres_   += repeat( eval( tvec_ ), R0 ) / scalar;
            trefres_ += repeat( eval( trefvec_ ), R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with addition assignment with evaluated vector (OP/s, compile time)
      {
         test_  = "Scaled repeat with addition assignment with evaluated vector (OP/s, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat<R0>( eval( vec_ ) ) / scalar;
            sres_   += repeat<R0>( eval( vec_ ) ) / scalar;
            refres_ += repeat<R0>( eval( refvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += repeat<R0>( eval( tvec_ ) ) / scalar;
            tsres_   += repeat<R0>( eval( tvec_ ) ) / scalar;
            trefres_ += repeat<R0>( eval( trefvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled repeat with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled repeat with subtraction assignment with the given vector (s*OP, runtime)
      {
         test_  = "Scaled repeat with subtraction assignment with the given vector (s*OP, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= scalar * repeat( vec_, R0 );
            sres_   -= scalar * repeat( vec_, R0 );
            refres_ -= scalar * repeat( refvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= scalar * repeat( tvec_, R0 );
            tsres_   -= scalar * repeat( tvec_, R0 );
            trefres_ -= scalar * repeat( trefvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with subtraction assignment with the given vector (s*OP, compile time)
      {
         test_  = "Scaled repeat with subtraction assignment with the given vector (s*OP, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= scalar * repeat<R0>( vec_ );
            sres_   -= scalar * repeat<R0>( vec_ );
            refres_ -= scalar * repeat<R0>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= scalar * repeat<R0>( tvec_ );
            tsres_   -= scalar * repeat<R0>( tvec_ );
            trefres_ -= scalar * repeat<R0>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with subtraction assignment with evaluated vector (s*OP, runtime)
      {
         test_  = "Scaled repeat with subtraction assignment with evaluated vector (s*OP, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= scalar * repeat( eval( vec_ ), R0 );
            sres_   -= scalar * repeat( eval( vec_ ), R0 );
            refres_ -= scalar * repeat( eval( refvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= scalar * repeat( eval( tvec_ ), R0 );
            tsres_   -= scalar * repeat( eval( tvec_ ), R0 );
            trefres_ -= scalar * repeat( eval( trefvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with subtraction assignment with evaluated vector (s*OP, compile time)
      {
         test_  = "Scaled repeat with subtraction assignment with evaluated vector (s*OP, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= scalar * repeat<R0>( eval( vec_ ) );
            sres_   -= scalar * repeat<R0>( eval( vec_ ) );
            refres_ -= scalar * repeat<R0>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= scalar * repeat<R0>( eval( tvec_ ) );
            tsres_   -= scalar * repeat<R0>( eval( tvec_ ) );
            trefres_ -= scalar * repeat<R0>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled repeat with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled repeat with subtraction assignment with the given vector (OP*s, runtime)
      {
         test_  = "Scaled repeat with subtraction assignment with the given vector (OP*s, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat( vec_, R0 ) * scalar;
            sres_   -= repeat( vec_, R0 ) * scalar;
            refres_ -= repeat( refvec_, R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= repeat( tvec_, R0 ) * scalar;
            tsres_   -= repeat( tvec_, R0 ) * scalar;
            trefres_ -= repeat( trefvec_, R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with subtraction assignment with the given vector (OP*s, compile time)
      {
         test_  = "Scaled repeat with subtraction assignment with the given vector (OP*s, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat<R0>( vec_ ) * scalar;
            sres_   -= repeat<R0>( vec_ ) * scalar;
            refres_ -= repeat<R0>( refvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= repeat<R0>( tvec_ ) * scalar;
            tsres_   -= repeat<R0>( tvec_ ) * scalar;
            trefres_ -= repeat<R0>( trefvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with subtraction assignment with evaluated vector (OP*s, runtime)
      {
         test_  = "Scaled repeat with subtraction assignment with evaluated vector (OP*s, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat( eval( vec_ ), R0 ) * scalar;
            sres_   -= repeat( eval( vec_ ), R0 ) * scalar;
            refres_ -= repeat( eval( refvec_ ), R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= repeat( eval( tvec_ ), R0 ) * scalar;
            tsres_   -= repeat( eval( tvec_ ), R0 ) * scalar;
            trefres_ -= repeat( eval( trefvec_ ), R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with subtraction assignment with evaluated vector (OP*s, compile time)
      {
         test_  = "Scaled repeat with subtraction assignment with evaluated vector (OP*s, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat<R0>( eval( vec_ ) ) * scalar;
            sres_   -= repeat<R0>( eval( vec_ ) ) * scalar;
            refres_ -= repeat<R0>( eval( refvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= repeat<R0>( eval( tvec_ ) ) * scalar;
            tsres_   -= repeat<R0>( eval( tvec_ ) ) * scalar;
            trefres_ -= repeat<R0>( eval( trefvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled repeat with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled repeat with subtraction assignment with the given vector (OP/s, runtime)
      {
         test_  = "Scaled repeat with subtraction assignment with the given vector (OP/s, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat( vec_, R0 ) / scalar;
            sres_   -= repeat( vec_, R0 ) / scalar;
            refres_ -= repeat( refvec_, R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= repeat( tvec_, R0 ) / scalar;
            tsres_   -= repeat( tvec_, R0 ) / scalar;
            trefres_ -= repeat( trefvec_, R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with subtraction assignment with the given vector (OP/s, compile time)
      {
         test_  = "Scaled repeat with subtraction assignment with the given vector (OP/s, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat<R0>( vec_ ) / scalar;
            sres_   -= repeat<R0>( vec_ ) / scalar;
            refres_ -= repeat<R0>( refvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= repeat<R0>( tvec_ ) / scalar;
            tsres_   -= repeat<R0>( tvec_ ) / scalar;
            trefres_ -= repeat<R0>( trefvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with subtraction assignment with evaluated vector (OP/s, runtime)
      {
         test_  = "Scaled repeat with subtraction assignment with evaluated vector (OP/s, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat( eval( vec_ ), R0 ) / scalar;
            sres_   -= repeat( eval( vec_ ), R0 ) / scalar;
            refres_ -= repeat( eval( refvec_ ), R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= repeat( eval( tvec_ ), R0 ) / scalar;
            tsres_   -= repeat( eval( tvec_ ), R0 ) / scalar;
            trefres_ -= repeat( eval( trefvec_ ), R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with subtraction assignment with evaluated vector (OP/s, compile time)
      {
         test_  = "Scaled repeat with subtraction assignment with evaluated vector (OP/s, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat<R0>( eval( vec_ ) ) / scalar;
            sres_   -= repeat<R0>( eval( vec_ ) ) / scalar;
            refres_ -= repeat<R0>( eval( refvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= repeat<R0>( eval( tvec_ ) ) / scalar;
            tsres_   -= repeat<R0>( eval( tvec_ ) ) / scalar;
            trefres_ -= repeat<R0>( eval( trefvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled repeat with multiplication assignment (s*OP)
      //=====================================================================================

      // Scaled repeat with multiplication assignment with the given vector (s*OP, runtime)
      {
         test_  = "Scaled repeat with multiplication assignment with the given vector (s*OP, runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= scalar * repeat( vec_, R0 );
            sres_   *= scalar * repeat( vec_, R0 );
            refres_ *= scalar * repeat( refvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= scalar * repeat( tvec_, R0 );
            tsres_   *= scalar * repeat( tvec_, R0 );
            trefres_ *= scalar * repeat( trefvec_, R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with multiplication assignment with the given vector (s*OP, compile time)
      {
         test_  = "Scaled repeat with multiplication assignment with the given vector (s*OP, compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= scalar * repeat<R0>( vec_ );
            sres_   *= scalar * repeat<R0>( vec_ );
            refres_ *= scalar * repeat<R0>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= scalar * repeat<R0>( tvec_ );
            tsres_   *= scalar * repeat<R0>( tvec_ );
            trefres_ *= scalar * repeat<R0>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with multiplication assignment with evaluated vector (s*OP, runtime)
      {
         test_  = "Scaled repeat with multiplication assignment with evaluated vector (s*OP, runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= scalar * repeat( eval( vec_ ), R0 );
            sres_   *= scalar * repeat( eval( vec_ ), R0 );
            refres_ *= scalar * repeat( eval( refvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= scalar * repeat( eval( tvec_ ), R0 );
            tsres_   *= scalar * repeat( eval( tvec_ ), R0 );
            trefres_ *= scalar * repeat( eval( trefvec_ ), R0 );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with multiplication assignment with evaluated vector (s*OP, compile time)
      {
         test_  = "Scaled repeat with multiplication assignment with evaluated vector (s*OP, compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= scalar * repeat<R0>( eval( vec_ ) );
            sres_   *= scalar * repeat<R0>( eval( vec_ ) );
            refres_ *= scalar * repeat<R0>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= scalar * repeat<R0>( eval( tvec_ ) );
            tsres_   *= scalar * repeat<R0>( eval( tvec_ ) );
            trefres_ *= scalar * repeat<R0>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled repeat with multiplication assignment (OP*s)
      //=====================================================================================

      // Scaled repeat with multiplication assignment with the given vector (OP*s, runtime)
      {
         test_  = "Scaled repeat with multiplication assignment with the given vector (OP*s, runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= repeat( vec_, R0 ) * scalar;
            sres_   *= repeat( vec_, R0 ) * scalar;
            refres_ *= repeat( refvec_, R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= repeat( tvec_, R0 ) * scalar;
            tsres_   *= repeat( tvec_, R0 ) * scalar;
            trefres_ *= repeat( trefvec_, R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with multiplication assignment with the given vector (OP*s, compile time)
      {
         test_  = "Scaled repeat with multiplication assignment with the given vector (OP*s, compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= repeat<R0>( vec_ ) * scalar;
            sres_   *= repeat<R0>( vec_ ) * scalar;
            refres_ *= repeat<R0>( refvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= repeat<R0>( tvec_ ) * scalar;
            tsres_   *= repeat<R0>( tvec_ ) * scalar;
            trefres_ *= repeat<R0>( trefvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with multiplication assignment with evaluated vector (OP*s, runtime)
      {
         test_  = "Scaled repeat with multiplication assignment with evaluated vector (OP*s, runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= repeat( eval( vec_ ), R0 ) * scalar;
            sres_   *= repeat( eval( vec_ ), R0 ) * scalar;
            refres_ *= repeat( eval( refvec_ ), R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= repeat( eval( tvec_ ), R0 ) * scalar;
            tsres_   *= repeat( eval( tvec_ ), R0 ) * scalar;
            trefres_ *= repeat( eval( trefvec_ ), R0 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with multiplication assignment with evaluated vector (OP*s, compile time)
      {
         test_  = "Scaled repeat with multiplication assignment with evaluated vector (OP*s, compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= repeat<R0>( eval( vec_ ) ) * scalar;
            sres_   *= repeat<R0>( eval( vec_ ) ) * scalar;
            refres_ *= repeat<R0>( eval( refvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= repeat<R0>( eval( tvec_ ) ) * scalar;
            tsres_   *= repeat<R0>( eval( tvec_ ) ) * scalar;
            trefres_ *= repeat<R0>( eval( trefvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled repeat with multiplication assignment (OP/s)
      //=====================================================================================

      // Scaled repeat with multiplication assignment with the given vector (OP/s, runtime)
      {
         test_  = "Scaled repeat with multiplication assignment with the given vector (OP/s, runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= repeat( vec_, R0 ) / scalar;
            sres_   *= repeat( vec_, R0 ) / scalar;
            refres_ *= repeat( refvec_, R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= repeat( tvec_, R0 ) / scalar;
            tsres_   *= repeat( tvec_, R0 ) / scalar;
            trefres_ *= repeat( trefvec_, R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with multiplication assignment with the given vector (OP/s, compile time)
      {
         test_  = "Scaled repeat with multiplication assignment with the given vector (OP/s, compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= repeat<R0>( vec_ ) / scalar;
            sres_   *= repeat<R0>( vec_ ) / scalar;
            refres_ *= repeat<R0>( refvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= repeat<R0>( tvec_ ) / scalar;
            tsres_   *= repeat<R0>( tvec_ ) / scalar;
            trefres_ *= repeat<R0>( trefvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with multiplication assignment with evaluated vector (OP/s, runtime)
      {
         test_  = "Scaled repeat with multiplication assignment with evaluated vector (OP/s, runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= repeat( eval( vec_ ), R0 ) / scalar;
            sres_   *= repeat( eval( vec_ ), R0 ) / scalar;
            refres_ *= repeat( eval( refvec_ ), R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= repeat( eval( tvec_ ), R0 ) / scalar;
            tsres_   *= repeat( eval( tvec_ ), R0 ) / scalar;
            trefres_ *= repeat( eval( trefvec_ ), R0 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled repeat with multiplication assignment with evaluated vector (OP/s, compile time)
      {
         test_  = "Scaled repeat with multiplication assignment with evaluated vector (OP/s, compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            dres_   *= repeat<R0>( eval( vec_ ) ) / scalar;
            sres_   *= repeat<R0>( eval( vec_ ) ) / scalar;
            refres_ *= repeat<R0>( eval( refvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   *= repeat<R0>( eval( tvec_ ) ) / scalar;
            tsres_   *= repeat<R0>( eval( tvec_ ) ) / scalar;
            trefres_ *= repeat<R0>( eval( trefvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the transpose dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the transpose vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANS_OPERATION > 1 )
   {
      using blaze::repeat;


      //=====================================================================================
      // Transpose repeat operation
      //=====================================================================================

      // Transpose repeat operation with the given vector (runtime)
      {
         test_  = "Transpose repeat operation with the given vector (runtime)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_   = trans( repeat( vec_, R0 ) );
            tsres_   = trans( repeat( vec_, R0 ) );
            trefres_ = trans( repeat( refvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = trans( repeat( tvec_, R0 ) );
            sres_   = trans( repeat( tvec_, R0 ) );
            refres_ = trans( repeat( trefvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose repeat operation with the given vector (compile time)
      {
         test_  = "Transpose repeat operation with the given vector (compile time)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_   = trans( repeat<R0>( vec_ ) );
            tsres_   = trans( repeat<R0>( vec_ ) );
            trefres_ = trans( repeat<R0>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = trans( repeat<R0>( tvec_ ) );
            sres_   = trans( repeat<R0>( tvec_ ) );
            refres_ = trans( repeat<R0>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose repeat operation with evaluated vector (runtime)
      {
         test_  = "Transpose repeat operation with evaluated vector (runtime)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_   = trans( repeat( eval( vec_ ), R0 ) );
            tsres_   = trans( repeat( eval( vec_ ), R0 ) );
            trefres_ = trans( repeat( eval( refvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = trans( repeat( eval( tvec_ ), R0 ) );
            sres_   = trans( repeat( eval( tvec_ ), R0 ) );
            refres_ = trans( repeat( eval( trefvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose repeat operation with evaluated vector (compile time)
      {
         test_  = "Transpose repeat operation with evaluated vector (compile time)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_   = trans( repeat<R0>( eval( vec_ ) ) );
            tsres_   = trans( repeat<R0>( eval( vec_ ) ) );
            trefres_ = trans( repeat<R0>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = trans( repeat<R0>( eval( tvec_ ) ) );
            sres_   = trans( repeat<R0>( eval( tvec_ ) ) );
            refres_ = trans( repeat<R0>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }


      //=====================================================================================
      // Transpose repeat with addition assignment
      //=====================================================================================

      // Transpose repeat with addition assignment with the given vector (runtime)
      {
         test_  = "Transpose repeat with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += trans( repeat( vec_, R0 ) );
            tsres_   += trans( repeat( vec_, R0 ) );
            trefres_ += trans( repeat( refvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += trans( repeat( tvec_, R0 ) );
            sres_   += trans( repeat( tvec_, R0 ) );
            refres_ += trans( repeat( trefvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose repeat with addition assignment with the given vector (compile time)
      {
         test_  = "Transpose repeat with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += trans( repeat<R0>( vec_ ) );
            tsres_   += trans( repeat<R0>( vec_ ) );
            trefres_ += trans( repeat<R0>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += trans( repeat<R0>( tvec_ ) );
            sres_   += trans( repeat<R0>( tvec_ ) );
            refres_ += trans( repeat<R0>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose repeat with addition assignment with evaluated vector (runtime)
      {
         test_  = "Transpose repeat with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += trans( repeat( eval( vec_ ), R0 ) );
            tsres_   += trans( repeat( eval( vec_ ), R0 ) );
            trefres_ += trans( repeat( eval( refvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += trans( repeat( eval( tvec_ ), R0 ) );
            sres_   += trans( repeat( eval( tvec_ ), R0 ) );
            refres_ += trans( repeat( eval( trefvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose repeat with addition assignment with evaluated vector (compile time)
      {
         test_  = "Transpose repeat with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += trans( repeat<R0>( eval( vec_ ) ) );
            tsres_   += trans( repeat<R0>( eval( vec_ ) ) );
            trefres_ += trans( repeat<R0>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += trans( repeat<R0>( eval( tvec_ ) ) );
            sres_   += trans( repeat<R0>( eval( tvec_ ) ) );
            refres_ += trans( repeat<R0>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }


      //=====================================================================================
      // Transpose repeat with subtraction assignment
      //=====================================================================================

      // Transpose repeat with subtraction assignment with the given vector (runtime)
      {
         test_  = "Transpose repeat with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= trans( repeat( vec_, R0 ) );
            tsres_   -= trans( repeat( vec_, R0 ) );
            trefres_ -= trans( repeat( refvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= trans( repeat( tvec_, R0 ) );
            sres_   -= trans( repeat( tvec_, R0 ) );
            refres_ -= trans( repeat( trefvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose repeat with subtraction assignment with the given vector (compile time)
      {
         test_  = "Transpose repeat with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= trans( repeat<R0>( vec_ ) );
            tsres_   -= trans( repeat<R0>( vec_ ) );
            trefres_ -= trans( repeat<R0>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= trans( repeat<R0>( tvec_ ) );
            sres_   -= trans( repeat<R0>( tvec_ ) );
            refres_ -= trans( repeat<R0>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose repeat with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Transpose repeat with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= trans( repeat( eval( vec_ ), R0 ) );
            tsres_   -= trans( repeat( eval( vec_ ), R0 ) );
            trefres_ -= trans( repeat( eval( refvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= trans( repeat( eval( tvec_ ), R0 ) );
            sres_   -= trans( repeat( eval( tvec_ ), R0 ) );
            refres_ -= trans( repeat( eval( trefvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose repeat with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Transpose repeat with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= trans( repeat<R0>( eval( vec_ ) ) );
            tsres_   -= trans( repeat<R0>( eval( vec_ ) ) );
            trefres_ -= trans( repeat<R0>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= trans( repeat<R0>( eval( tvec_ ) ) );
            sres_   -= trans( repeat<R0>( eval( tvec_ ) ) );
            refres_ -= trans( repeat<R0>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }


      //=====================================================================================
      // Transpose repeat with multiplication assignment
      //=====================================================================================

      // Transpose repeat with multiplication assignment with the given vector (runtime)
      {
         test_  = "Transpose repeat with multiplication assignment with the given vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initTransposeResults();
            tdres_   *= trans( repeat( vec_, R0 ) );
            tsres_   *= trans( repeat( vec_, R0 ) );
            trefres_ *= trans( repeat( refvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   *= trans( repeat( tvec_, R0 ) );
            sres_   *= trans( repeat( tvec_, R0 ) );
            refres_ *= trans( repeat( trefvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose repeat with multiplication assignment with the given vector (compile time)
      {
         test_  = "Transpose repeat with multiplication assignment with the given vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initTransposeResults();
            tdres_   *= trans( repeat<R0>( vec_ ) );
            tsres_   *= trans( repeat<R0>( vec_ ) );
            trefres_ *= trans( repeat<R0>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   *= trans( repeat<R0>( tvec_ ) );
            sres_   *= trans( repeat<R0>( tvec_ ) );
            refres_ *= trans( repeat<R0>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose repeat with multiplication assignment with evaluated vector (runtime)
      {
         test_  = "Transpose repeat with multiplication assignment with evaluated vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initTransposeResults();
            tdres_   *= trans( repeat( eval( vec_ ), R0 ) );
            tsres_   *= trans( repeat( eval( vec_ ), R0 ) );
            trefres_ *= trans( repeat( eval( refvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   *= trans( repeat( eval( tvec_ ), R0 ) );
            sres_   *= trans( repeat( eval( tvec_ ), R0 ) );
            refres_ *= trans( repeat( eval( trefvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose repeat with multiplication assignment with evaluated vector (compile time)
      {
         test_  = "Transpose repeat with multiplication assignment with evaluated vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initTransposeResults();
            tdres_   *= trans( repeat<R0>( eval( vec_ ) ) );
            tsres_   *= trans( repeat<R0>( eval( vec_ ) ) );
            trefres_ *= trans( repeat<R0>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   *= trans( repeat<R0>( eval( tvec_ ) ) );
            sres_   *= trans( repeat<R0>( eval( tvec_ ) ) );
            refres_ *= trans( repeat<R0>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the conjugate transpose dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the conjugate transpose vector repeat operation with plain assignment,
// addition assignment, subtraction assignment, multiplication assignment, and division assignment.
// In case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testCTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION > 1 )
   {
      using blaze::repeat;


      //=====================================================================================
      // Conjugate transpose repeat operation
      //=====================================================================================

      // Conjugate transpose repeat operation with the given vector (runtime)
      {
         test_  = "Conjugate transpose repeat operation with the given vector (runtime)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( repeat( vec_, R0 ) );
            tsres_   = ctrans( repeat( vec_, R0 ) );
            trefres_ = ctrans( repeat( refvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = ctrans( repeat( tvec_, R0 ) );
            sres_   = ctrans( repeat( tvec_, R0 ) );
            refres_ = ctrans( repeat( trefvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose repeat operation with the given vector (compile time)
      {
         test_  = "Conjugate transpose repeat operation with the given vector (compile time)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( repeat<R0>( vec_ ) );
            tsres_   = ctrans( repeat<R0>( vec_ ) );
            trefres_ = ctrans( repeat<R0>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = ctrans( repeat<R0>( tvec_ ) );
            sres_   = ctrans( repeat<R0>( tvec_ ) );
            refres_ = ctrans( repeat<R0>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose repeat operation with evaluated vector (runtime)
      {
         test_  = "Conjugate transpose repeat operation with evaluated vector (runtime)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( repeat( eval( vec_ ), R0 ) );
            tsres_   = ctrans( repeat( eval( vec_ ), R0 ) );
            trefres_ = ctrans( repeat( eval( refvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = ctrans( repeat( eval( tvec_ ), R0 ) );
            sres_   = ctrans( repeat( eval( tvec_ ), R0 ) );
            refres_ = ctrans( repeat( eval( trefvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose repeat operation with evaluated vector (compile time)
      {
         test_  = "Conjugate transpose repeat operation with evaluated vector (compile time)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( repeat<R0>( eval( vec_ ) ) );
            tsres_   = ctrans( repeat<R0>( eval( vec_ ) ) );
            trefres_ = ctrans( repeat<R0>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = ctrans( repeat<R0>( eval( tvec_ ) ) );
            sres_   = ctrans( repeat<R0>( eval( tvec_ ) ) );
            refres_ = ctrans( repeat<R0>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }


      //=====================================================================================
      // Conjugate transpose repeat with addition assignment
      //=====================================================================================

      // Conjugate transpose repeat with addition assignment with the given vector (runtime)
      {
         test_  = "Conjugate transpose repeat with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += ctrans( repeat( vec_, R0 ) );
            tsres_   += ctrans( repeat( vec_, R0 ) );
            trefres_ += ctrans( repeat( refvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += ctrans( repeat( tvec_, R0 ) );
            sres_   += ctrans( repeat( tvec_, R0 ) );
            refres_ += ctrans( repeat( trefvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose repeat with addition assignment with the given vector (compile time)
      {
         test_  = "Conjugate transpose repeat with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += ctrans( repeat<R0>( vec_ ) );
            tsres_   += ctrans( repeat<R0>( vec_ ) );
            trefres_ += ctrans( repeat<R0>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += ctrans( repeat<R0>( tvec_ ) );
            sres_   += ctrans( repeat<R0>( tvec_ ) );
            refres_ += ctrans( repeat<R0>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose repeat with addition assignment with evaluated vector (runtime)
      {
         test_  = "Conjugate transpose repeat with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += ctrans( repeat( eval( vec_ ), R0 ) );
            tsres_   += ctrans( repeat( eval( vec_ ), R0 ) );
            trefres_ += ctrans( repeat( eval( refvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += ctrans( repeat( eval( tvec_ ), R0 ) );
            sres_   += ctrans( repeat( eval( tvec_ ), R0 ) );
            refres_ += ctrans( repeat( eval( trefvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose repeat with addition assignment with evaluated vector (compile time)
      {
         test_  = "Conjugate transpose repeat with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += ctrans( repeat<R0>( eval( vec_ ) ) );
            tsres_   += ctrans( repeat<R0>( eval( vec_ ) ) );
            trefres_ += ctrans( repeat<R0>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += ctrans( repeat<R0>( eval( tvec_ ) ) );
            sres_   += ctrans( repeat<R0>( eval( tvec_ ) ) );
            refres_ += ctrans( repeat<R0>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }


      //=====================================================================================
      // Conjugate transpose repeat with subtraction assignment
      //=====================================================================================

      // Conjugate transpose repeat with subtraction assignment with the given vector (runtime)
      {
         test_  = "Conjugate transpose repeat with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= ctrans( repeat( vec_, R0 ) );
            tsres_   -= ctrans( repeat( vec_, R0 ) );
            trefres_ -= ctrans( repeat( refvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= ctrans( repeat( tvec_, R0 ) );
            sres_   -= ctrans( repeat( tvec_, R0 ) );
            refres_ -= ctrans( repeat( trefvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose repeat with subtraction assignment with the given vector (compile time)
      {
         test_  = "Conjugate transpose repeat with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= ctrans( repeat<R0>( vec_ ) );
            tsres_   -= ctrans( repeat<R0>( vec_ ) );
            trefres_ -= ctrans( repeat<R0>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= ctrans( repeat<R0>( tvec_ ) );
            sres_   -= ctrans( repeat<R0>( tvec_ ) );
            refres_ -= ctrans( repeat<R0>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose repeat with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Conjugate transpose repeat with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= ctrans( repeat( eval( vec_ ), R0 ) );
            tsres_   -= ctrans( repeat( eval( vec_ ), R0 ) );
            trefres_ -= ctrans( repeat( eval( refvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= ctrans( repeat( eval( tvec_ ), R0 ) );
            sres_   -= ctrans( repeat( eval( tvec_ ), R0 ) );
            refres_ -= ctrans( repeat( eval( trefvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose repeat with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Conjugate transpose repeat with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= ctrans( repeat<R0>( eval( vec_ ) ) );
            tsres_   -= ctrans( repeat<R0>( eval( vec_ ) ) );
            trefres_ -= ctrans( repeat<R0>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= ctrans( repeat<R0>( eval( tvec_ ) ) );
            sres_   -= ctrans( repeat<R0>( eval( tvec_ ) ) );
            refres_ -= ctrans( repeat<R0>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }


      //=====================================================================================
      // Conjugate transpose repeat with multiplication assignment
      //=====================================================================================

      // Conjugate transpose repeat with multiplication assignment with the given vector (runtime)
      {
         test_  = "Conjugate transpose repeat with multiplication assignment with the given vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initTransposeResults();
            tdres_   *= ctrans( repeat( vec_, R0 ) );
            tsres_   *= ctrans( repeat( vec_, R0 ) );
            trefres_ *= ctrans( repeat( refvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   *= ctrans( repeat( tvec_, R0 ) );
            sres_   *= ctrans( repeat( tvec_, R0 ) );
            refres_ *= ctrans( repeat( trefvec_, R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose repeat with multiplication assignment with the given vector (compile time)
      {
         test_  = "Conjugate transpose repeat with multiplication assignment with the given vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initTransposeResults();
            tdres_   *= ctrans( repeat<R0>( vec_ ) );
            tsres_   *= ctrans( repeat<R0>( vec_ ) );
            trefres_ *= ctrans( repeat<R0>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   *= ctrans( repeat<R0>( tvec_ ) );
            sres_   *= ctrans( repeat<R0>( tvec_ ) );
            refres_ *= ctrans( repeat<R0>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose repeat with multiplication assignment with evaluated vector (runtime)
      {
         test_  = "Conjugate transpose repeat with multiplication assignment with evaluated vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initTransposeResults();
            tdres_   *= ctrans( repeat( eval( vec_ ), R0 ) );
            tsres_   *= ctrans( repeat( eval( vec_ ), R0 ) );
            trefres_ *= ctrans( repeat( eval( refvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   *= ctrans( repeat( eval( tvec_ ), R0 ) );
            sres_   *= ctrans( repeat( eval( tvec_ ), R0 ) );
            refres_ *= ctrans( repeat( eval( trefvec_ ), R0 ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose repeat with multiplication assignment with evaluated vector (compile time)
      {
         test_  = "Conjugate transpose repeat with multiplication assignment with evaluated vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initTransposeResults();
            tdres_   *= ctrans( repeat<R0>( eval( vec_ ) ) );
            tsres_   *= ctrans( repeat<R0>( eval( vec_ ) ) );
            trefres_ *= ctrans( repeat<R0>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   *= ctrans( repeat<R0>( eval( tvec_ ) ) );
            sres_   *= ctrans( repeat<R0>( eval( tvec_ ) ) );
            refres_ *= ctrans( repeat<R0>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the abs dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the abs vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testAbsOperation()
{
#if BLAZETEST_MATHTEST_TEST_ABS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ABS_OPERATION > 1 )
   {
      testCustomOperation( blaze::Abs(), "abs" );
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the conjugate dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the conjugate vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testConjOperation()
{
#if BLAZETEST_MATHTEST_TEST_CONJ_OPERATION
   if( BLAZETEST_MATHTEST_TEST_CONJ_OPERATION > 1 )
   {
      testCustomOperation( blaze::Conj(), "conj" );
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the \a real dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the \a real vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testRealOperation()
{
#if BLAZETEST_MATHTEST_TEST_REAL_OPERATION
   if( BLAZETEST_MATHTEST_TEST_REAL_OPERATION > 1 )
   {
      testCustomOperation( blaze::Real(), "real" );
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the \a imag dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the \a imag vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testImagOperation()
{
#if BLAZETEST_MATHTEST_TEST_IMAG_OPERATION
   if( BLAZETEST_MATHTEST_TEST_IMAG_OPERATION > 1 )
   {
      testCustomOperation( blaze::Imag(), "imag" );
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the evaluated dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the evaluated vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testEvalOperation()
{
#if BLAZETEST_MATHTEST_TEST_EVAL_OPERATION
   if( BLAZETEST_MATHTEST_TEST_EVAL_OPERATION > 1 )
   {
      testCustomOperation( blaze::Eval(), "eval" );
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the serialized dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the serialized vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testSerialOperation()
{
#if BLAZETEST_MATHTEST_TEST_SERIAL_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SERIAL_OPERATION > 1 )
   {
      testCustomOperation( blaze::Serial(), "serial" );
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the non-aliased dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the non-aliased vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testNoAliasOperation()
{
#if BLAZETEST_MATHTEST_TEST_NOALIAS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NOALIAS_OPERATION > 1 )
   {
      testCustomOperation( blaze::NoAlias(), "noalias" );
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the non-SIMD dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the non-SIMD vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testNoSIMDOperation()
{
#if BLAZETEST_MATHTEST_TEST_NOSIMD_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NOSIMD_OPERATION > 1 )
   {
      testCustomOperation( blaze::NoSIMD(), "nosimd" );
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the subvector-wise dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the subvector-wise vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testSubvectorOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION > 1 )
   {
      using blaze::repeat;


      if( vec_.size() == 0UL )
         return;


      //=====================================================================================
      // Subvector-wise repeat operation
      //=====================================================================================

      // Subvector-wise repeat operation with the given vector (runtime)
      {
         test_  = "Subvector-wise repeat operation with the given vector (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) = subvector( repeat( vec_   , R0 ), index, size );
               subvector( sres_  , index, size ) = subvector( repeat( vec_   , R0 ), index, size );
               subvector( refres_, index, size ) = subvector( repeat( refvec_, R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  , index, size ) = subvector( repeat( tvec_   , R0 ), index, size );
               subvector( tsres_  , index, size ) = subvector( repeat( tvec_   , R0 ), index, size );
               subvector( trefres_, index, size ) = subvector( repeat( trefvec_, R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Subvector-wise repeat operation with the given vector (compile time)
      {
         test_  = "Subvector-wise repeat operation with the given vector (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) = subvector( repeat<R0>( vec_    ), index, size );
               subvector( sres_  , index, size ) = subvector( repeat<R0>( vec_    ), index, size );
               subvector( refres_, index, size ) = subvector( repeat<R0>( refvec_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  ,  index, size ) = subvector( repeat<R0>( tvec_    ), index, size );
               subvector( tsres_  ,  index, size ) = subvector( repeat<R0>( tvec_    ), index, size );
               subvector( trefres_,  index, size ) = subvector( repeat<R0>( trefvec_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Subvector-wise repeat operation with evaluated vector (runtime)
      {
         test_  = "Subvector-wise repeat operation with evaluated vector (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) = subvector( repeat( eval( vec_    ), R0 ), index, size );
               subvector( sres_  , index, size ) = subvector( repeat( eval( vec_    ), R0 ), index, size );
               subvector( refres_, index, size ) = subvector( repeat( eval( refvec_ ), R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  , index, size ) = subvector( repeat( eval( tvec_    ), R0 ), index, size );
               subvector( tsres_  , index, size ) = subvector( repeat( eval( tvec_    ), R0 ), index, size );
               subvector( trefres_, index, size ) = subvector( repeat( eval( trefvec_ ), R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Subvector-wise repeat operation with evaluated vector (compile time)
      {
         test_  = "Subvector-wise repeat operation with evaluated vector (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) = subvector( repeat<R0>( eval( vec_    ) ), index, size );
               subvector( sres_  , index, size ) = subvector( repeat<R0>( eval( vec_    ) ), index, size );
               subvector( refres_, index, size ) = subvector( repeat<R0>( eval( refvec_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  , index, size ) = subvector( repeat<R0>( eval( tvec_    ) ), index, size );
               subvector( tsres_  , index, size ) = subvector( repeat<R0>( eval( tvec_    ) ), index, size );
               subvector( trefres_, index, size ) = subvector( repeat<R0>( eval( trefvec_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Subvector-wise repeat with addition assignment
      //=====================================================================================

      // Subvector-wise repeat with addition assignment with the given vector (runtime)
      {
         test_  = "Subvector-wise repeat with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) += subvector( repeat( vec_   , R0 ), index, size );
               subvector( sres_  , index, size ) += subvector( repeat( vec_   , R0 ), index, size );
               subvector( refres_, index, size ) += subvector( repeat( refvec_, R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  , index, size ) += subvector( repeat( tvec_   , R0 ), index, size );
               subvector( tsres_  , index, size ) += subvector( repeat( tvec_   , R0 ), index, size );
               subvector( trefres_, index, size ) += subvector( repeat( trefvec_, R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Subvector-wise repeat with addition assignment with the given vector (compile time)
      {
         test_  = "Subvector-wise repeat with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) += subvector( repeat<R0>( vec_    ), index, size );
               subvector( sres_  , index, size ) += subvector( repeat<R0>( vec_    ), index, size );
               subvector( refres_, index, size ) += subvector( repeat<R0>( refvec_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  ,  index, size ) += subvector( repeat<R0>( tvec_    ), index, size );
               subvector( tsres_  ,  index, size ) += subvector( repeat<R0>( tvec_    ), index, size );
               subvector( trefres_,  index, size ) += subvector( repeat<R0>( trefvec_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Subvector-wise repeat with addition assignment with evaluated vector (runtime)
      {
         test_  = "Subvector-wise repeat with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) += subvector( repeat( eval( vec_    ), R0 ), index, size );
               subvector( sres_  , index, size ) += subvector( repeat( eval( vec_    ), R0 ), index, size );
               subvector( refres_, index, size ) += subvector( repeat( eval( refvec_ ), R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  , index, size ) += subvector( repeat( eval( tvec_    ), R0 ), index, size );
               subvector( tsres_  , index, size ) += subvector( repeat( eval( tvec_    ), R0 ), index, size );
               subvector( trefres_, index, size ) += subvector( repeat( eval( trefvec_ ), R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Subvector-wise repeat with addition assignment with evaluated vector (compile time)
      {
         test_  = "Subvector-wise repeat with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) += subvector( repeat<R0>( eval( vec_    ) ), index, size );
               subvector( sres_  , index, size ) += subvector( repeat<R0>( eval( vec_    ) ), index, size );
               subvector( refres_, index, size ) += subvector( repeat<R0>( eval( refvec_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  , index, size ) += subvector( repeat<R0>( eval( tvec_    ) ), index, size );
               subvector( tsres_  , index, size ) += subvector( repeat<R0>( eval( tvec_    ) ), index, size );
               subvector( trefres_, index, size ) += subvector( repeat<R0>( eval( trefvec_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Subvector-wise repeat with subtraction assignment
      //=====================================================================================

      // Subvector-wise repeat with subtraction assignment with the given vector (runtime)
      {
         test_  = "Subvector-wise repeat with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) -= subvector( repeat( vec_   , R0 ), index, size );
               subvector( sres_  , index, size ) -= subvector( repeat( vec_   , R0 ), index, size );
               subvector( refres_, index, size ) -= subvector( repeat( refvec_, R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  , index, size ) -= subvector( repeat( tvec_   , R0 ), index, size );
               subvector( tsres_  , index, size ) -= subvector( repeat( tvec_   , R0 ), index, size );
               subvector( trefres_, index, size ) -= subvector( repeat( trefvec_, R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Subvector-wise repeat with subtraction assignment with the given vector (compile time)
      {
         test_  = "Subvector-wise repeat with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) -= subvector( repeat<R0>( vec_    ), index, size );
               subvector( sres_  , index, size ) -= subvector( repeat<R0>( vec_    ), index, size );
               subvector( refres_, index, size ) -= subvector( repeat<R0>( refvec_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  ,  index, size ) -= subvector( repeat<R0>( tvec_    ), index, size );
               subvector( tsres_  ,  index, size ) -= subvector( repeat<R0>( tvec_    ), index, size );
               subvector( trefres_,  index, size ) -= subvector( repeat<R0>( trefvec_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Subvector-wise repeat with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Subvector-wise repeat with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) -= subvector( repeat( eval( vec_    ), R0 ), index, size );
               subvector( sres_  , index, size ) -= subvector( repeat( eval( vec_    ), R0 ), index, size );
               subvector( refres_, index, size ) -= subvector( repeat( eval( refvec_ ), R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  , index, size ) -= subvector( repeat( eval( tvec_    ), R0 ), index, size );
               subvector( tsres_  , index, size ) -= subvector( repeat( eval( tvec_    ), R0 ), index, size );
               subvector( trefres_, index, size ) -= subvector( repeat( eval( trefvec_ ), R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Subvector-wise repeat with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Subvector-wise repeat with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) -= subvector( repeat<R0>( eval( vec_    ) ), index, size );
               subvector( sres_  , index, size ) -= subvector( repeat<R0>( eval( vec_    ) ), index, size );
               subvector( refres_, index, size ) -= subvector( repeat<R0>( eval( refvec_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  , index, size ) -= subvector( repeat<R0>( eval( tvec_    ) ), index, size );
               subvector( tsres_  , index, size ) -= subvector( repeat<R0>( eval( tvec_    ) ), index, size );
               subvector( trefres_, index, size ) -= subvector( repeat<R0>( eval( trefvec_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Subvector-wise repeat with multiplication assignment
      //=====================================================================================

      // Subvector-wise repeat with multiplication assignment with the given vector (runtime)
      {
         test_  = "Subvector-wise repeat with multiplication assignment with the given vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) *= subvector( repeat( vec_   , R0 ), index, size );
               subvector( sres_  , index, size ) *= subvector( repeat( vec_   , R0 ), index, size );
               subvector( refres_, index, size ) *= subvector( repeat( refvec_, R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  , index, size ) *= subvector( repeat( tvec_   , R0 ), index, size );
               subvector( tsres_  , index, size ) *= subvector( repeat( tvec_   , R0 ), index, size );
               subvector( trefres_, index, size ) *= subvector( repeat( trefvec_, R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Subvector-wise repeat with multiplication assignment with the given vector (compile time)
      {
         test_  = "Subvector-wise repeat with multiplication assignment with the given vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) *= subvector( repeat<R0>( vec_    ), index, size );
               subvector( sres_  , index, size ) *= subvector( repeat<R0>( vec_    ), index, size );
               subvector( refres_, index, size ) *= subvector( repeat<R0>( refvec_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  ,  index, size ) *= subvector( repeat<R0>( tvec_    ), index, size );
               subvector( tsres_  ,  index, size ) *= subvector( repeat<R0>( tvec_    ), index, size );
               subvector( trefres_,  index, size ) *= subvector( repeat<R0>( trefvec_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Subvector-wise repeat with multiplication assignment with evaluated vector (runtime)
      {
         test_  = "Subvector-wise repeat with multiplication assignment with evaluated vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) *= subvector( repeat( eval( vec_    ), R0 ), index, size );
               subvector( sres_  , index, size ) *= subvector( repeat( eval( vec_    ), R0 ), index, size );
               subvector( refres_, index, size ) *= subvector( repeat( eval( refvec_ ), R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  , index, size ) *= subvector( repeat( eval( tvec_    ), R0 ), index, size );
               subvector( tsres_  , index, size ) *= subvector( repeat( eval( tvec_    ), R0 ), index, size );
               subvector( trefres_, index, size ) *= subvector( repeat( eval( trefvec_ ), R0 ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Subvector-wise repeat with multiplication assignment with evaluated vector (compile time)
      {
         test_  = "Subvector-wise repeat with multiplication assignment with evaluated vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<dres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, dres_.size() - index );
               subvector( dres_  , index, size ) *= subvector( repeat<R0>( eval( vec_    ) ), index, size );
               subvector( sres_  , index, size ) *= subvector( repeat<R0>( eval( vec_    ) ), index, size );
               subvector( refres_, index, size ) *= subvector( repeat<R0>( eval( refvec_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tdres_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tdres_.size() - index );
               subvector( tdres_  , index, size ) *= subvector( repeat<R0>( eval( tvec_    ) ), index, size );
               subvector( tsres_  , index, size ) *= subvector( repeat<R0>( eval( tvec_    ) ), index, size );
               subvector( trefres_, index, size ) *= subvector( repeat<R0>( eval( trefvec_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      try {
         auto sv = subvector( repeat( vec_, R0 ), 1UL, vec_.size()*R0 );

         std::ostringstream oss;
         oss << " Test: Subvector construction\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense vector type:\n"
             << "     " << typeid( VT ).name() << "\n"
             << "   Result:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ex )
      {
         if( std::strcmp( ex.what(), "Invalid subvector specification" ) != 0 ) {
            std::ostringstream oss;
            oss << " Test: Subvector construction\n"
                << " Error: Wrong error message\n"
                << " Details:\n"
                << "   Error message: \"" << ex.what() << "\"\n"
                << "   Expected error message: \"Invalid subvector specification\"\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         auto sv = subvector( repeat( vec_, R0 ), vec_.size()*R0, 1UL );

         std::ostringstream oss;
         oss << " Test: Subvector construction\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense vector type:\n"
             << "     " << typeid( VT ).name() << "\n"
             << "   Result:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ex )
      {
         if( std::strcmp( ex.what(), "Invalid subvector specification" ) != 0 ) {
            std::ostringstream oss;
            oss << " Test: Subvector construction\n"
                << " Error: Wrong error message\n"
                << " Details:\n"
                << "   Error message: \"" << ex.what() << "\"\n"
                << "   Expected error message: \"Invalid subvector specification\"\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the subvector-wise dense vector repeat operation.
//
// \return void
//
// This function is called in case the subvector-wise vector repeat operation is not available
// for the given vector type \a VT.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testSubvectorOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the elements-wise dense vector repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the elements-wise vector repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment.
// In case any error resulting from the repeat operation or the subsequent assignment is
// detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testElementsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_ELEMENTS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ELEMENTS_OPERATION > 1 )
   {
      using blaze::repeat;


      if( vec_.size() == 0UL )
         return;


      std::vector<size_t> indices( vec_.size() * R0 );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );


      //=====================================================================================
      // Elements-wise repeat operation
      //=====================================================================================

      // Elements-wise repeat operation with the given vector (runtime)
      {
         test_  = "Elements-wise repeat operation with the given vector (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) = elements( repeat( vec_   , R0 ), &indices[index], n );
               elements( sres_  , &indices[index], n ) = elements( repeat( vec_   , R0 ), &indices[index], n );
               elements( refres_, &indices[index], n ) = elements( repeat( refvec_, R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) = elements( repeat( tvec_   , R0 ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) = elements( repeat( tvec_   , R0 ), &indices[index], n );
               elements( trefres_, &indices[index], n ) = elements( repeat( trefvec_, R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Elements-wise repeat operation with the given vector (compile time)
      {
         test_  = "Elements-wise repeat operation with the given vector (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) = elements( repeat<R0>( vec_    ), &indices[index], n );
               elements( sres_  , &indices[index], n ) = elements( repeat<R0>( vec_    ), &indices[index], n );
               elements( refres_, &indices[index], n ) = elements( repeat<R0>( refvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  ,  &indices[index], n ) = elements( repeat<R0>( tvec_    ), &indices[index], n );
               elements( tsres_  ,  &indices[index], n ) = elements( repeat<R0>( tvec_    ), &indices[index], n );
               elements( trefres_,  &indices[index], n ) = elements( repeat<R0>( trefvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Elements-wise repeat operation with evaluated vector (runtime)
      {
         test_  = "Elements-wise repeat operation with evaluated vector (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) = elements( repeat( eval( vec_    ), R0 ), &indices[index], n );
               elements( sres_  , &indices[index], n ) = elements( repeat( eval( vec_    ), R0 ), &indices[index], n );
               elements( refres_, &indices[index], n ) = elements( repeat( eval( refvec_ ), R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) = elements( repeat( eval( tvec_    ), R0 ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) = elements( repeat( eval( tvec_    ), R0 ), &indices[index], n );
               elements( trefres_, &indices[index], n ) = elements( repeat( eval( trefvec_ ), R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Elements-wise repeat operation with evaluated vector (compile time)
      {
         test_  = "Elements-wise repeat operation with evaluated vector (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) = elements( repeat<R0>( eval( vec_    ) ), &indices[index], n );
               elements( sres_  , &indices[index], n ) = elements( repeat<R0>( eval( vec_    ) ), &indices[index], n );
               elements( refres_, &indices[index], n ) = elements( repeat<R0>( eval( refvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) = elements( repeat<R0>( eval( tvec_    ) ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) = elements( repeat<R0>( eval( tvec_    ) ), &indices[index], n );
               elements( trefres_, &indices[index], n ) = elements( repeat<R0>( eval( trefvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Elements-wise repeat with addition assignment
      //=====================================================================================

      // Elements-wise repeat with addition assignment with the given vector (runtime)
      {
         test_  = "Elements-wise repeat with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) += elements( repeat( vec_   , R0 ), &indices[index], n );
               elements( sres_  , &indices[index], n ) += elements( repeat( vec_   , R0 ), &indices[index], n );
               elements( refres_, &indices[index], n ) += elements( repeat( refvec_, R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) += elements( repeat( tvec_   , R0 ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) += elements( repeat( tvec_   , R0 ), &indices[index], n );
               elements( trefres_, &indices[index], n ) += elements( repeat( trefvec_, R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Elements-wise repeat with addition assignment with the given vector (compile time)
      {
         test_  = "Elements-wise repeat with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) += elements( repeat<R0>( vec_    ), &indices[index], n );
               elements( sres_  , &indices[index], n ) += elements( repeat<R0>( vec_    ), &indices[index], n );
               elements( refres_, &indices[index], n ) += elements( repeat<R0>( refvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  ,  &indices[index], n ) += elements( repeat<R0>( tvec_    ), &indices[index], n );
               elements( tsres_  ,  &indices[index], n ) += elements( repeat<R0>( tvec_    ), &indices[index], n );
               elements( trefres_,  &indices[index], n ) += elements( repeat<R0>( trefvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Elements-wise repeat with addition assignment with evaluated vector (runtime)
      {
         test_  = "Elements-wise repeat with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) += elements( repeat( eval( vec_    ), R0 ), &indices[index], n );
               elements( sres_  , &indices[index], n ) += elements( repeat( eval( vec_    ), R0 ), &indices[index], n );
               elements( refres_, &indices[index], n ) += elements( repeat( eval( refvec_ ), R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) += elements( repeat( eval( tvec_    ), R0 ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) += elements( repeat( eval( tvec_    ), R0 ), &indices[index], n );
               elements( trefres_, &indices[index], n ) += elements( repeat( eval( trefvec_ ), R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Elements-wise repeat with addition assignment with evaluated vector (compile time)
      {
         test_  = "Elements-wise repeat with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) += elements( repeat<R0>( eval( vec_    ) ), &indices[index], n );
               elements( sres_  , &indices[index], n ) += elements( repeat<R0>( eval( vec_    ) ), &indices[index], n );
               elements( refres_, &indices[index], n ) += elements( repeat<R0>( eval( refvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) += elements( repeat<R0>( eval( tvec_    ) ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) += elements( repeat<R0>( eval( tvec_    ) ), &indices[index], n );
               elements( trefres_, &indices[index], n ) += elements( repeat<R0>( eval( trefvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Elements-wise repeat with subtraction assignment
      //=====================================================================================

      // Elements-wise repeat with subtraction assignment with the given vector (runtime)
      {
         test_  = "Elements-wise repeat with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) -= elements( repeat( vec_   , R0 ), &indices[index], n );
               elements( sres_  , &indices[index], n ) -= elements( repeat( vec_   , R0 ), &indices[index], n );
               elements( refres_, &indices[index], n ) -= elements( repeat( refvec_, R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) -= elements( repeat( tvec_   , R0 ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) -= elements( repeat( tvec_   , R0 ), &indices[index], n );
               elements( trefres_, &indices[index], n ) -= elements( repeat( trefvec_, R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Elements-wise repeat with subtraction assignment with the given vector (compile time)
      {
         test_  = "Elements-wise repeat with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) -= elements( repeat<R0>( vec_    ), &indices[index], n );
               elements( sres_  , &indices[index], n ) -= elements( repeat<R0>( vec_    ), &indices[index], n );
               elements( refres_, &indices[index], n ) -= elements( repeat<R0>( refvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  ,  &indices[index], n ) -= elements( repeat<R0>( tvec_    ), &indices[index], n );
               elements( tsres_  ,  &indices[index], n ) -= elements( repeat<R0>( tvec_    ), &indices[index], n );
               elements( trefres_,  &indices[index], n ) -= elements( repeat<R0>( trefvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Elements-wise repeat with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Elements-wise repeat with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) -= elements( repeat( eval( vec_    ), R0 ), &indices[index], n );
               elements( sres_  , &indices[index], n ) -= elements( repeat( eval( vec_    ), R0 ), &indices[index], n );
               elements( refres_, &indices[index], n ) -= elements( repeat( eval( refvec_ ), R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) -= elements( repeat( eval( tvec_    ), R0 ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) -= elements( repeat( eval( tvec_    ), R0 ), &indices[index], n );
               elements( trefres_, &indices[index], n ) -= elements( repeat( eval( trefvec_ ), R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Elements-wise repeat with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Elements-wise repeat with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) -= elements( repeat<R0>( eval( vec_    ) ), &indices[index], n );
               elements( sres_  , &indices[index], n ) -= elements( repeat<R0>( eval( vec_    ) ), &indices[index], n );
               elements( refres_, &indices[index], n ) -= elements( repeat<R0>( eval( refvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) -= elements( repeat<R0>( eval( tvec_    ) ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) -= elements( repeat<R0>( eval( tvec_    ) ), &indices[index], n );
               elements( trefres_, &indices[index], n ) -= elements( repeat<R0>( eval( trefvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Elements-wise repeat with multiplication assignment
      //=====================================================================================

      // Elements-wise repeat with multiplication assignment with the given vector (runtime)
      {
         test_  = "Elements-wise repeat with multiplication assignment with the given vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) *= elements( repeat( vec_   , R0 ), &indices[index], n );
               elements( sres_  , &indices[index], n ) *= elements( repeat( vec_   , R0 ), &indices[index], n );
               elements( refres_, &indices[index], n ) *= elements( repeat( refvec_, R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) *= elements( repeat( tvec_   , R0 ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) *= elements( repeat( tvec_   , R0 ), &indices[index], n );
               elements( trefres_, &indices[index], n ) *= elements( repeat( trefvec_, R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Elements-wise repeat with multiplication assignment with the given vector (compile time)
      {
         test_  = "Elements-wise repeat with multiplication assignment with the given vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) *= elements( repeat<R0>( vec_    ), &indices[index], n );
               elements( sres_  , &indices[index], n ) *= elements( repeat<R0>( vec_    ), &indices[index], n );
               elements( refres_, &indices[index], n ) *= elements( repeat<R0>( refvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  ,  &indices[index], n ) *= elements( repeat<R0>( tvec_    ), &indices[index], n );
               elements( tsres_  ,  &indices[index], n ) *= elements( repeat<R0>( tvec_    ), &indices[index], n );
               elements( trefres_,  &indices[index], n ) *= elements( repeat<R0>( trefvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Elements-wise repeat with multiplication assignment with evaluated vector (runtime)
      {
         test_  = "Elements-wise repeat with multiplication assignment with evaluated vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) *= elements( repeat( eval( vec_    ), R0 ), &indices[index], n );
               elements( sres_  , &indices[index], n ) *= elements( repeat( eval( vec_    ), R0 ), &indices[index], n );
               elements( refres_, &indices[index], n ) *= elements( repeat( eval( refvec_ ), R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) *= elements( repeat( eval( tvec_    ), R0 ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) *= elements( repeat( eval( tvec_    ), R0 ), &indices[index], n );
               elements( trefres_, &indices[index], n ) *= elements( repeat( eval( trefvec_ ), R0 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Elements-wise repeat with multiplication assignment with evaluated vector (compile time)
      {
         test_  = "Elements-wise repeat with multiplication assignment with evaluated vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) *= elements( repeat<R0>( eval( vec_    ) ), &indices[index], n );
               elements( sres_  , &indices[index], n ) *= elements( repeat<R0>( eval( vec_    ) ), &indices[index], n );
               elements( refres_, &indices[index], n ) *= elements( repeat<R0>( eval( refvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) *= elements( repeat<R0>( eval( tvec_    ) ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) *= elements( repeat<R0>( eval( tvec_    ) ), &indices[index], n );
               elements( trefres_, &indices[index], n ) *= elements( repeat<R0>( eval( trefvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      try {
         auto e =  elements( repeat( vec_, R0 ), blaze::index_sequence<128*R0>() );

         std::ostringstream oss;
         oss << " Test: Elements construction (index_sequence)\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense vector type:\n"
             << "     " << typeid( VT ).name() << "\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ex )
      {
         if( std::strcmp( ex.what(), "Invalid element access index" ) != 0 ) {
            std::ostringstream oss;
            oss << " Test: Elements construction (index_sequence)\n"
                << " Error: Wrong error message\n"
                << " Details:\n"
                << "   Error message: \"" << ex.what() << "\"\n"
                << "   Expected error message: \"Invalid element access index\"\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         auto e = elements( repeat( vec_, R0 ), { vec_.size()*R0 } );

         std::ostringstream oss;
         oss << " Test: Elements construction (initializer_list)\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense vector type:\n"
             << "     " << typeid( VT ).name() << "\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ex )
      {
         if( std::strcmp( ex.what(), "Invalid element access index" ) != 0 ) {
            std::ostringstream oss;
            oss << " Test: Elements construction (initializer_list)\n"
                << " Error: Wrong error message\n"
                << " Details:\n"
                << "   Error message: \"" << ex.what() << "\"\n"
                << "   Expected error message: \"Invalid element access index\"\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         auto e = elements( repeat( vec_, R0 ), [index=vec_.size()*R0]( size_t ){ return index; }, 1UL );

         std::ostringstream oss;
         oss << " Test: Elements construction (lambda)\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense vector type:\n"
             << "     " << typeid( VT ).name() << "\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ex )
      {
         if( std::strcmp( ex.what(), "Invalid element access index" ) != 0 ) {
            std::ostringstream oss;
            oss << " Test: Elements construction (lambda)\n"
                << " Error: Wrong error message\n"
                << " Details:\n"
                << "   Error message: \"" << ex.what() << "\"\n"
                << "   Expected error message: \"Invalid element access index\"\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the elements-wise dense vector repeat operation.
//
// \return void
//
// This function is called in case the elements-wise vector repeat operation is not available
// for the given vector type \a VT.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::testElementsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the customized dense vector repeat operation.
//
// \param op The custom operation to be tested.
// \param name The human-readable name of the operation.
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the vector repeat operation with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment in combination with
// a custom operation. In case any error resulting from the repeat operation or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
template< typename OP >   // Type of the custom operation
void OperationTest<VT,R0>::testCustomOperation( OP op, const std::string& name )
{
   using blaze::repeat;


   //=====================================================================================
   // Customized repeat operation
   //=====================================================================================

   // Customized repeat operation with the given vector (runtime)
   {
      test_  = "Customized repeat operation with the given vector (runtime)";
      error_ = "Failed repeat operation";

      try {
         initResults();
         dres_   = op( repeat( vec_, R0 ) );
         sres_   = op( repeat( vec_, R0 ) );
         refres_ = op( repeat( refvec_, R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   = op( repeat( tvec_, R0 ) );
         tsres_   = op( repeat( tvec_, R0 ) );
         trefres_ = op( repeat( trefvec_, R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized repeat operation with the given vector (compile time)
   {
      test_  = "Customized repeat operation with the given vector (compile time)";
      error_ = "Failed repeat operation";

      try {
         initResults();
         dres_   = op( repeat<R0>( vec_ ) );
         sres_   = op( repeat<R0>( vec_ ) );
         refres_ = op( repeat<R0>( refvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   = op( repeat<R0>( tvec_ ) );
         tsres_   = op( repeat<R0>( tvec_ ) );
         trefres_ = op( repeat<R0>( trefvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized repeat operation with evaluated vector (runtime)
   {
      test_  = "Customized repeat operation with evaluated vector (runtime)";
      error_ = "Failed repeat operation";

      try {
         initResults();
         dres_   = op( repeat( eval( vec_ ), R0 ) );
         sres_   = op( repeat( eval( vec_ ), R0 ) );
         refres_ = op( repeat( eval( refvec_ ), R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   = op( repeat( eval( tvec_ ), R0 ) );
         tsres_   = op( repeat( eval( tvec_ ), R0 ) );
         trefres_ = op( repeat( eval( trefvec_ ), R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized repeat operation with evaluated vector (compile time)
   {
      test_  = "Customized repeat operation with evaluated vector (compile time)";
      error_ = "Failed repeat operation";

      try {
         initResults();
         dres_   = op( repeat<R0>( eval( vec_ ) ) );
         sres_   = op( repeat<R0>( eval( vec_ ) ) );
         refres_ = op( repeat<R0>( eval( refvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   = op( repeat<R0>( eval( tvec_ ) ) );
         tsres_   = op( repeat<R0>( eval( tvec_ ) ) );
         trefres_ = op( repeat<R0>( eval( trefvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }


   //=====================================================================================
   // Customized repeat with addition assignment
   //=====================================================================================

   // Customized repeat with addition assignment with the given vector (runtime)
   {
      test_  = "Customized repeat with addition assignment with the given vector (runtime)";
      error_ = "Failed addition assignment";

      try {
         initResults();
         dres_   += op( repeat( vec_, R0 ) );
         sres_   += op( repeat( vec_, R0 ) );
         refres_ += op( repeat( refvec_, R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   += op( repeat( tvec_, R0 ) );
         tsres_   += op( repeat( tvec_, R0 ) );
         trefres_ += op( repeat( trefvec_, R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized repeat with addition assignment with the given vector (compile time)
   {
      test_  = "Customized repeat with addition assignment with the given vector (compile time)";
      error_ = "Failed addition assignment";

      try {
         initResults();
         dres_   += op( repeat<R0>( vec_ ) );
         sres_   += op( repeat<R0>( vec_ ) );
         refres_ += op( repeat<R0>( refvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   += op( repeat<R0>( tvec_ ) );
         tsres_   += op( repeat<R0>( tvec_ ) );
         trefres_ += op( repeat<R0>( trefvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized repeat with addition assignment with evaluated vector (runtime)
   {
      test_  = "Customized repeat with addition assignment with evaluated vector (runtime)";
      error_ = "Failed addition assignment";

      try {
         initResults();
         dres_   += op( repeat( eval( vec_ ), R0 ) );
         sres_   += op( repeat( eval( vec_ ), R0 ) );
         refres_ += op( repeat( eval( refvec_ ), R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   += op( repeat( eval( tvec_ ), R0 ) );
         tsres_   += op( repeat( eval( tvec_ ), R0 ) );
         trefres_ += op( repeat( eval( trefvec_ ), R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized repeat with addition assignment with evaluated vector (compile time)
   {
      test_  = "Customized repeat with addition assignment with evaluated vector (compile time)";
      error_ = "Failed addition assignment";

      try {
         initResults();
         dres_   += op( repeat<R0>( eval( vec_ ) ) );
         sres_   += op( repeat<R0>( eval( vec_ ) ) );
         refres_ += op( repeat<R0>( eval( refvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   += op( repeat<R0>( eval( tvec_ ) ) );
         tsres_   += op( repeat<R0>( eval( tvec_ ) ) );
         trefres_ += op( repeat<R0>( eval( trefvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }


   //=====================================================================================
   // Customized repeat with subtraction assignment
   //=====================================================================================

   // Customized repeat with subtraction assignment with the given vector (runtime)
   {
      test_  = "Customized repeat with subtraction assignment with the given vector (runtime)";
      error_ = "Failed subtraction assignment";

      try {
         initResults();
         dres_   -= op( repeat( vec_, R0 ) );
         sres_   -= op( repeat( vec_, R0 ) );
         refres_ -= op( repeat( refvec_, R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   -= op( repeat( tvec_, R0 ) );
         tsres_   -= op( repeat( tvec_, R0 ) );
         trefres_ -= op( repeat( trefvec_, R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized repeat with subtraction assignment with the given vector (compile time)
   {
      test_  = "Customized repeat with subtraction assignment with the given vector (compile time)";
      error_ = "Failed subtraction assignment";

      try {
         initResults();
         dres_   -= op( repeat<R0>( vec_ ) );
         sres_   -= op( repeat<R0>( vec_ ) );
         refres_ -= op( repeat<R0>( refvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   -= op( repeat<R0>( tvec_ ) );
         tsres_   -= op( repeat<R0>( tvec_ ) );
         trefres_ -= op( repeat<R0>( trefvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized repeat with subtraction assignment with evaluated vector (runtime)
   {
      test_  = "Customized repeat with subtraction assignment with evaluated vector (runtime)";
      error_ = "Failed subtraction assignment";

      try {
         initResults();
         dres_   -= op( repeat( eval( vec_ ), R0 ) );
         sres_   -= op( repeat( eval( vec_ ), R0 ) );
         refres_ -= op( repeat( eval( refvec_ ), R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   -= op( repeat( eval( tvec_ ), R0 ) );
         tsres_   -= op( repeat( eval( tvec_ ), R0 ) );
         trefres_ -= op( repeat( eval( trefvec_ ), R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized repeat with subtraction assignment with evaluated vector (compile time)
   {
      test_  = "Customized repeat with subtraction assignment with evaluated vector (compile time)";
      error_ = "Failed subtraction assignment";

      try {
         initResults();
         dres_   -= op( repeat<R0>( eval( vec_ ) ) );
         sres_   -= op( repeat<R0>( eval( vec_ ) ) );
         refres_ -= op( repeat<R0>( eval( refvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   -= op( repeat<R0>( eval( tvec_ ) ) );
         tsres_   -= op( repeat<R0>( eval( tvec_ ) ) );
         trefres_ -= op( repeat<R0>( eval( trefvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }


   //=====================================================================================
   // Customized repeat with multiplication assignment
   //=====================================================================================

   // Customized repeat with multiplication assignment with the given vector (runtime)
   {
      test_  = "Customized repeat with multiplication assignment with the given vector (runtime)";
      error_ = "Failed multiplication assignment";

      try {
         initResults();
         dres_   *= op( repeat( vec_, R0 ) );
         sres_   *= op( repeat( vec_, R0 ) );
         refres_ *= op( repeat( refvec_, R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   *= op( repeat( tvec_, R0 ) );
         tsres_   *= op( repeat( tvec_, R0 ) );
         trefres_ *= op( repeat( trefvec_, R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized repeat with multiplication assignment with the given vector (compile time)
   {
      test_  = "Customized repeat with multiplication assignment with the given vector (compile time)";
      error_ = "Failed multiplication assignment";

      try {
         initResults();
         dres_   *= op( repeat<R0>( vec_ ) );
         sres_   *= op( repeat<R0>( vec_ ) );
         refres_ *= op( repeat<R0>( refvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   *= op( repeat<R0>( tvec_ ) );
         tsres_   *= op( repeat<R0>( tvec_ ) );
         trefres_ *= op( repeat<R0>( trefvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized repeat with multiplication assignment with evaluated vector (runtime)
   {
      test_  = "Customized repeat with multiplication assignment with evaluated vector (runtime)";
      error_ = "Failed multiplication assignment";

      try {
         initResults();
         dres_   *= op( repeat( eval( vec_ ), R0 ) );
         sres_   *= op( repeat( eval( vec_ ), R0 ) );
         refres_ *= op( repeat( eval( refvec_ ), R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   *= op( repeat( eval( tvec_ ), R0 ) );
         tsres_   *= op( repeat( eval( tvec_ ), R0 ) );
         trefres_ *= op( repeat( eval( trefvec_ ), R0 ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized repeat with multiplication assignment with evaluated vector (compile time)
   {
      test_  = "Customized repeat with multiplication assignment with evaluated vector (compile time)";
      error_ = "Failed multiplication assignment";

      try {
         initResults();
         dres_   *= op( repeat<R0>( eval( vec_ ) ) );
         sres_   *= op( repeat<R0>( eval( vec_ ) ) );
         refres_ *= op( repeat<R0>( eval( refvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   *= op( repeat<R0>( eval( tvec_ ) ) );
         tsres_   *= op( repeat<R0>( eval( tvec_ ) ) );
         trefres_ *= op( repeat<R0>( eval( trefvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }
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
// \exception std::runtime_error Incorrect dense result detected.
// \exception std::runtime_error Incorrect sparse result detected.
//
// This function is called after each test case to check and compare the computed results. The
// two template arguments \a LT and \a RT indicate the types of the left-hand side and right-hand
// side operands used for the computations.
*/
template< typename VT   // Type of the dense vector
        , size_t R0 >   // Compile time repetitions
template< typename T >  // Type of the vector operand
void OperationTest<VT,R0>::checkResults()
{
   using blaze::IsRowVector;

   if( !isEqual( dres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result vector detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense " << ( IsRowVector<T>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( T ).name() << "\n"
          << "   Result:\n" << dres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( sres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result vector detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense " << ( IsRowVector<T>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( T ).name() << "\n"
          << "   Result:\n" << sres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking and comparing the computed transpose results.
//
// \return void
// \exception std::runtime_error Incorrect dense result detected.
// \exception std::runtime_error Incorrect sparse result detected.
//
// This function is called after each test case to check and compare the computed transpose
// results. The two template arguments \a LT and \a RT indicate the types of the left-hand
// side and right-hand side operands used for the computations.
*/
template< typename VT   // Type of the dense vector
        , size_t R0 >   // Compile time repetitions
template< typename T >  // Type of the vector operand
void OperationTest<VT,R0>::checkTransposeResults()
{
   using blaze::IsRowVector;

   if( !isEqual( tdres_, trefres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result vector detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense " << ( IsRowVector<T>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( T ).name() << "\n"
          << "   Transpose result:\n" << tdres_ << "\n"
          << "   Expected result:\n" << trefres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( tsres_, trefres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result vector detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense " << ( IsRowVector<T>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( T ).name() << "\n"
          << "   Transpose result:\n" << tsres_ << "\n"
          << "   Expected result:\n" << trefres_ << "\n";
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
/*!\brief Initializing the non-transpose result vectors.
//
// \return void
//
// This function is called before each non-transpose test case to initialize the according result
// vectors to random values.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::initResults()
{
   const blaze::UnderlyingBuiltin_t<DRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<DRE> max( randmax );

   resize( dres_, size( vec_ )*R0 );
   randomize( dres_, min, max );

   sres_   = dres_;
   refres_ = dres_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initializing the transpose result vectors.
//
// \return void
//
// This function is called before each transpose test case to initialize the according result
// vectors to random values.
*/
template< typename VT  // Type of the dense vector
        , size_t R0 >  // Compile time repetitions
void OperationTest<VT,R0>::initTransposeResults()
{
   const blaze::UnderlyingBuiltin_t<TDRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<TDRE> max( randmax );

   resize( tdres_, size( tvec_ )*R0 );
   randomize( tdres_, min, max );

   tsres_   = tdres_;
   trefres_ = tdres_;
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
// test. The two template arguments \a LT and \a RT indicate the types of the left-hand side and
// right-hand side operands used for the computations.
*/
template< typename VT   // Type of the dense vector
        , size_t R0 >   // Compile time repetitions
template< typename T >  // Type of the vector operand
void OperationTest<VT,R0>::convertException( const std::exception& ex )
{
   using blaze::IsRowVector;

   std::ostringstream oss;
   oss << " Test : " << test_ << "\n"
       << " Error: " << error_ << "\n"
       << " Details:\n"
       << "   Random seed = " << blaze::getSeed() << "\n"
       << "   Dense " << ( IsRowVector<T>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
       << "     " << typeid( T ).name() << "\n"
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
/*!\brief Testing the repeat operation for a specific vector type.
//
// \param creator The creator for the dense vector.
// \return void
*/
template< typename VT >  // Type of the dense vector
void runTest( const Creator<VT>& creator )
{
   for( size_t rep=0UL; rep<BLAZETEST_REPETITIONS; ++rep ) {
      OperationTest<VT,3UL>{ creator };
      OperationTest<VT,6UL>{ creator };
      OperationTest<VT,7UL>{ creator };
      OperationTest<VT,16UL>{ creator };
      OperationTest<VT,17UL>{ creator };
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the definition of a dense vector repeat operation test case.
*/
#define DEFINE_DVECREPEAT_OPERATION_TEST( VT ) \
   extern template class blazetest::mathtest::operations::dvecrepeat::OperationTest<VT>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a dense vector repeat operation test case.
*/
#define RUN_DVECREPEAT_OPERATION_TEST( C ) \
   blazetest::mathtest::operations::dvecrepeat::runTest( C )
/*! \endcond */
//*************************************************************************************************

} // namespace dvecrepeat

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
