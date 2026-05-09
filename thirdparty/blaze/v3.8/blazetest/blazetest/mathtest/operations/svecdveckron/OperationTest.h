//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/svecdveckron/OperationTest.h
//  \brief Header file for the sparse vector/dense vector Kronecker product operation test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_SVECDVECKRON_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_SVECDVECKRON_OPERATIONTEST_H_


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
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/Functors.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/traits/KronTrait.h>
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

namespace svecdveckron {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the sparse vector/dense vector Kronecker product operation test.
//
// This class template represents one particular vector Kronecker product test between two vectors
// of a particular type. The two template arguments \a VT1 and \a VT2 represent the types of the
// left-hand side and right-hand side vector, respectively.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
class OperationTest
{
 private:
   //**Enumerations********************************************************************************
   enum { TF = blaze::IsRowVector<VT1>::value };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   using ET1 = blaze::ElementType_t<VT1>;  //!< Element type 1
   using ET2 = blaze::ElementType_t<VT2>;  //!< Element type 2

   using TVT1 = blaze::TransposeType_t<VT1>;  //!< Transpose vector type 1
   using TVT2 = blaze::TransposeType_t<VT2>;  //!< Transpose vector type 2

   using SRE  = blaze::KronTrait_t<VT1,VT2>;    //!< Sparse result type
   using TSRE = blaze::KronTrait_t<TVT1,TVT2>;  //!< Transpose sparse result type
   using SET  = blaze::ElementType_t<SRE>;      //!< Element type of the sparse result

   using DRE  = blaze::DynamicVector<SET,TF>;  //!< Dense result type
   using TDRE = blaze::TransposeType_t<DRE>;   //!< Transpose dense result type
   using DET  = blaze::ElementType_t<DRE>;     //!< Element type of the dense result

   using RT1 = blaze::DynamicVector<ET1,TF>;  //!< Reference type 1
   using RT2 = blaze::DynamicVector<ET2,TF>;  //!< Reference type 2
   using RRE = blaze::KronTrait_t<RT1,RT2>;   //!< Reference result type

   using TRT1 = blaze::TransposeType_t<RT1>;    //!< Transpose reference type 1
   using TRT2 = blaze::TransposeType_t<RT2>;    //!< Transpose reference type 2
   using TRRE = blaze::KronTrait_t<TRT1,TRT2>;  //!< Transpose reference result type

   //! Type of the vector/vector Kronecker product expression
   using VecVecKronExprType =
      blaze::RemoveCVRef_t< decltype( kron( std::declval<VT1>(), std::declval<VT2>() ) ) >;

   //! Type of the transpose vector/transpose vector Kronecker product expression
   using TVecTVecKronExprType =
      blaze::RemoveCVRef_t< decltype( kron( std::declval<TVT1>(), std::declval<TVT2>() ) ) >;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit OperationTest( const Creator<VT1>& creator1, const Creator<VT2>& creator2 );
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
   template< typename LT, typename RT > void checkResults();
   template< typename LT, typename RT > void checkTransposeResults();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initResults();
   void initTransposeResults();
   template< typename LT, typename RT > void convertException( const std::exception& ex );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VT1  lhs_;      //!< The left-hand side sparse vector.
   VT2  rhs_;      //!< The right-hand side dense vector.
   DRE  dres_;     //!< The dense vector for the result of the vector Kronecker product.
   SRE  sres_;     //!< The sparse vector for the result of the vector Kronecker product.
   RT1  reflhs_;   //!< The reference left-hand side vector.
   RT2  refrhs_;   //!< The reference right-hand side vector.
   RRE  refres_;   //!< The reference result.
   TVT1 tlhs_;     //!< The transpose left-hand side vector.
   TVT2 trhs_;     //!< The transpose right-hand side vector.
   TDRE tdres_;    //!< The dense vector for the result of the transpose vector Kronecker product.
   TSRE tsres_;    //!< The sparse vector for the result of the transpose vector Kronecker product.
   TRT1 treflhs_;  //!< The reference left-hand side transpose vector.
   TRT2 trefrhs_;  //!< The reference right-hand side transpose vector.
   TRRE trefres_;  //!< The transpose reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT1  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT2  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TVT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TVT2 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( RT1  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( RT2  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( RRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TRT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TRT2 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TRRE );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( DRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( SRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TDRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TSRE );

   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , VT2  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , RT1  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TVT1, TVT2 );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TVT1, TRT1 );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( RT1 , RT2  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TRT1, TRT2 );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , RRE  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , DRE  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , SRE  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TVT1, TRRE );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TVT1, TDRE );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TVT1, TSRE );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET1, blaze::ElementType_t<TVT1>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET2, blaze::ElementType_t<TVT2>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<DRE>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<TDRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<SRE>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<SRE>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<TSRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<DRE>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT1, blaze::TransposeType_t<TVT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT2, blaze::TransposeType_t<TVT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RT1, blaze::TransposeType_t<TRT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RT2, blaze::TransposeType_t<TRT2> );

   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG     ( VecVecKronExprType, blaze::ResultType_t<VecVecKronExprType>    );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_DIFFERENT_TRANSPOSE_FLAG( VecVecKronExprType, blaze::TransposeType_t<VecVecKronExprType> );

   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG     ( TVecTVecKronExprType, blaze::ResultType_t<TVecTVecKronExprType>    );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_DIFFERENT_TRANSPOSE_FLAG( TVecTVecKronExprType, blaze::TransposeType_t<TVecTVecKronExprType> );

   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( VecVecKronExprType  , blaze::BaseType_t<VecVecKronExprType  > );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TVecTVecKronExprType, blaze::BaseType_t<TVecTVecKronExprType> );
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
/*!\brief Constructor for the sparse vector/dense vector Kronecker product operation test.
//
// \param creator1 The creator for the left-hand side sparse vector of the vector Kronecker product.
// \param creator2 The creator for the right-hand side dense vector of the vector Kronecker product.
// \exception std::runtime_error Operation error detected.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
OperationTest<VT1,VT2>::OperationTest( const Creator<VT1>& creator1, const Creator<VT2>& creator2 )
   : lhs_( creator1() )    // The left-hand side sparse vector
   , rhs_( creator2() )    // The right-hand side dense vector
   , dres_()               // The dense vector for the result of the vector Kronecker product
   , sres_()               // The sparse vector for the result of the vector Kronecker product
   , reflhs_( lhs_ )       // The reference left-hand side vector
   , refrhs_( rhs_ )       // The reference right-hand side vector
   , refres_()             // The reference result
   , tlhs_( trans(lhs_) )  // The transpose left-hand side vector
   , trhs_( trans(rhs_) )  // The transpose right-hand side vector
   , tdres_()              // The dense vector for the result of the transpose vector Kronecker product
   , tsres_()              // The sparse vector for the result of the transpose vector Kronecker product
   , treflhs_( tlhs_ )     // The reference left-hand side transpose vector
   , trefrhs_( trhs_ )     // The reference right-hand side transpose vector
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
/*!\brief Tests on the initial status of the vectors.
//
// \return void
// \exception std::runtime_error Initialization error detected.
//
// This function runs tests on the initial status of the vectors. In case any initialization
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testInitialStatus()
{
   //=====================================================================================
   // Performing initial tests with the given vectors
   //=====================================================================================

   // Checking the size of the left-hand side operand
   if( lhs_.size() != reflhs_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side sparse operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Detected size = " << lhs_.size() << "\n"
          << "   Expected size = " << reflhs_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the size of the right-hand side operand
   if( rhs_.size() != refrhs_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side dense operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
          << "   Detected size = " << rhs_.size() << "\n"
          << "   Expected size = " << refrhs_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the left-hand side operand
   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of left-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Current initialization:\n" << lhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the right-hand side operand
   if( !isEqual( rhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of right-hand side dense operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing initial tests with the transpose types
   //=====================================================================================

   // Checking the size of the left-hand side operand
   if( tlhs_.size() != treflhs_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of transpose left-hand side sparse operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose sparse vector type:\n"
          << "     " << typeid( TVT1 ).name() << "\n"
          << "   Detected size = " << tlhs_.size() << "\n"
          << "   Expected size = " << treflhs_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the size of the right-hand side operand
   if( trhs_.size() != trefrhs_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of transpose right-hand side dense operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Detected size = " << trhs_.size() << "\n"
          << "   Expected size = " << trefrhs_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the left-hand side operand
   if( !isEqual( tlhs_, treflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of transpose left-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose sparse vector type:\n"
          << "     " << typeid( TVT1 ).name() << "\n"
          << "   Current initialization:\n" << tlhs_ << "\n"
          << "   Expected initialization:\n" << treflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the right-hand side operand
   if( !isEqual( trhs_, trefrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of transpose right-hand side dense operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Current initialization:\n" << trhs_ << "\n"
          << "   Expected initialization:\n" << trefrhs_ << "\n";
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testAssignment()
{
   //=====================================================================================
   // Performing an assignment with the given vectors
   //=====================================================================================

   try {
      lhs_ = reflhs_;
      rhs_ = refrhs_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the given vectors\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side dense vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of left-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Current initialization:\n" << lhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( rhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of right-hand side dense operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing an assignment with the transpose types
   //=====================================================================================

   try {
      tlhs_ = treflhs_;
      trhs_ = trefrhs_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the transpose types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose left-hand side sparse vector type:\n"
          << "     " << typeid( TVT1 ).name() << "\n"
          << "   Transpose right-hand side dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( tlhs_, treflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of transpose left-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose sparse vector type:\n"
          << "     " << typeid( TVT1 ).name() << "\n"
          << "   Current initialization:\n" << tlhs_ << "\n"
          << "   Expected initialization:\n" << treflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( trhs_, trefrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of transpose right-hand side dense operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Current initialization:\n" << trhs_ << "\n"
          << "   Expected initialization:\n" << trefrhs_ << "\n";
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testEvaluation()
{
   using blaze::IsRowVector;


   //=====================================================================================
   // Testing the evaluation with the given vectors
   //=====================================================================================

   {
      const auto res   ( evaluate( kron( lhs_   , rhs_    ) ) );
      const auto refres( evaluate( kron( reflhs_, refrhs_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given vectors\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse " << ( IsRowVector<VT1>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side dense " << ( IsRowVector<VT2>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
             << "     " << typeid( rhs_ ).name() << "\n"
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
      const auto res   ( evaluate( kron( eval( lhs_ )   , eval( rhs_ )    ) ) );
      const auto refres( evaluate( kron( eval( reflhs_ ), eval( refrhs_ ) ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated vectors\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse " << ( IsRowVector<VT1>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side dense " << ( IsRowVector<VT2>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
             << "     " << typeid( rhs_ ).name() << "\n"
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
   // Testing the evaluation with the transpose types
   //=====================================================================================

   {
      const auto res   ( evaluate( kron( tlhs_   , trhs_    ) ) );
      const auto refres( evaluate( kron( treflhs_, trefrhs_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the transpose vectors\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse " << ( IsRowVector<TVT1>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
             << "     " << typeid( tlhs_ ).name() << "\n"
             << "   Right-hand side dense " << ( IsRowVector<TVT2>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
             << "     " << typeid( trhs_ ).name() << "\n"
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
      const auto res   ( evaluate( kron( eval( tlhs_ )   , eval( trhs_ )    ) ) );
      const auto refres( evaluate( kron( eval( treflhs_ ), eval( trefrhs_ ) ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated transpose vectors\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse " << ( IsRowVector<TVT1>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
             << "     " << typeid( tlhs_ ).name() << "\n"
             << "   Right-hand side dense " << ( IsRowVector<TVT2>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
             << "     " << typeid( trhs_ ).name() << "\n"
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testElementAccess()
{
   using blaze::equal;


   //=====================================================================================
   // Testing the element access with the given vectors
   //=====================================================================================

   if( lhs_.size() > 0UL && rhs_.size() > 0UL )
   {
      const size_t n( lhs_.size() - 1UL );

      if( !equal( kron( lhs_, rhs_ )[n], kron( reflhs_, refrhs_ )[n] ) ||
          !equal( kron( lhs_, rhs_ ).at(n), kron( reflhs_, refrhs_ ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of Kronecker product expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( kron( lhs_, eval( rhs_ ) )[n], kron( reflhs_, eval( refrhs_ ) )[n] ) ||
          !equal( kron( lhs_, eval( rhs_ ) ).at(n), kron( reflhs_, eval( refrhs_ ) ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated Kronecker product expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( kron( eval( lhs_ ), rhs_ )[n], kron( eval( reflhs_ ), refrhs_ )[n] ) ||
          !equal( kron( eval( lhs_ ), rhs_ ).at(n), kron( eval( reflhs_ ), refrhs_ ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated Kronecker product expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( kron( eval( lhs_ ), eval( rhs_ ) )[n], kron( eval( reflhs_ ), eval( refrhs_ ) )[n] ) |
          !equal( kron( eval( lhs_ ), eval( rhs_ ) ).at(n), kron( eval( reflhs_ ), eval( refrhs_ ) ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated Kronecker product expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      kron( lhs_, rhs_ ).at( lhs_.size() );

      std::ostringstream oss;
      oss << " Test : Checked element access of Kronecker product expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side sparse vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side dense vector type:\n"
          << "     " << typeid( VT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}


   //=====================================================================================
   // Testing the element access with the transpose types
   //=====================================================================================

   if( tlhs_.size() > 0UL && trhs_.size() > 0UL )
   {
      const size_t n( tlhs_.size() - 1UL );

      if( !equal( kron( tlhs_, trhs_ )[n], kron( treflhs_, trefrhs_ )[n] ) ||
          !equal( kron( tlhs_, trhs_ ).at(n), kron( treflhs_, trefrhs_ ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of transpose Kronecker product expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Transpose left-hand side sparse vector type:\n"
             << "     " << typeid( TVT1 ).name() << "\n"
             << "   Transpose right-hand side dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( kron( tlhs_, eval( trhs_ ) )[n], kron( treflhs_, eval( trefrhs_ ) )[n] ) ||
          !equal( kron( tlhs_, eval( trhs_ ) ).at(n), kron( treflhs_, eval( trefrhs_ ) ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated transpose Kronecker product expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Transpose left-hand side sparse vector type:\n"
             << "     " << typeid( TVT1 ).name() << "\n"
             << "   Transpose right-hand side dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( kron( eval( tlhs_ ), trhs_ )[n], kron( eval( treflhs_ ), trefrhs_ )[n] ) ||
          !equal( kron( eval( tlhs_ ), trhs_ ).at(n), kron( eval( treflhs_ ), trefrhs_ ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated transpose Kronecker product expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Transpose left-hand side sparse vector type:\n"
             << "     " << typeid( TVT1 ).name() << "\n"
             << "   Transpose right-hand side dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( kron( eval( tlhs_ ), eval( trhs_ ) )[n], kron( eval( treflhs_ ), eval( trefrhs_ ) )[n] ) ||
          !equal( kron( eval( tlhs_ ), eval( trhs_ ) ).at(n), kron( eval( treflhs_ ), eval( trefrhs_ ) ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated transpose Kronecker product expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Transpose left-hand side sparse vector type:\n"
             << "     " << typeid( TVT1 ).name() << "\n"
             << "   Transpose right-hand side dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      kron( tlhs_, trhs_ ).at( tlhs_.size() );

      std::ostringstream oss;
      oss << " Test : Checked element access of transpose Kronecker product expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose left-hand side sparse vector type:\n"
          << "     " << typeid( TVT1 ).name() << "\n"
          << "   Transpose right-hand side dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the plain vector Kronecker product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the Kronecker product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Kronecker product
      //=====================================================================================

      // Kronecker product with the given vectors
      {
         test_  = "Kronecker product with the given vectors";
         error_ = "Failed Kronecker product operation";

         try {
            initResults();
            dres_   = kron( lhs_, rhs_ );
            sres_   = kron( lhs_, rhs_ );
            refres_ = kron( reflhs_, refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = kron( tlhs_, trhs_ );
            tsres_   = kron( tlhs_, trhs_ );
            trefres_ = kron( treflhs_, trefrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Kronecker product with evaluated vectors
      {
         test_  = "Kronecker product with evaluated vectors";
         error_ = "Failed Kronecker product operation";

         try {
            initResults();
            dres_   = kron( eval( lhs_ ), eval( rhs_ ) );
            sres_   = kron( eval( lhs_ ), eval( rhs_ ) );
            refres_ = kron( eval( reflhs_ ), eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = kron( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   = kron( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ = kron( eval( treflhs_ ), eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Kronecker product with addition assignment
      //=====================================================================================

      // Kronecker product with addition assignment with the given vectors
      {
         test_  = "Kronecker product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += kron( lhs_, rhs_ );
            sres_   += kron( lhs_, rhs_ );
            refres_ += kron( reflhs_, refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += kron( tlhs_, trhs_ );
            tsres_   += kron( tlhs_, trhs_ );
            trefres_ += kron( treflhs_, trefrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Kronecker product with addition assignment with the given vectors
      {
         test_  = "Kronecker product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += kron( eval( lhs_ ), eval( rhs_ ) );
            sres_   += kron( eval( lhs_ ), eval( rhs_ ) );
            refres_ += kron( eval( reflhs_ ), eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += kron( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   += kron( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ += kron( eval( treflhs_ ), eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Kronecker product with subtraction assignment
      //=====================================================================================

      // Kronecker product with subtraction assignment with the given vectors
      {
         test_  = "Kronecker product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= kron( lhs_, rhs_ );
            sres_   -= kron( lhs_, rhs_ );
            refres_ -= kron( reflhs_, refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= kron( tlhs_, trhs_ );
            tsres_   -= kron( tlhs_, trhs_ );
            trefres_ -= kron( treflhs_, trefrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Kronecker product with subtraction assignment with evaluated vectors
      {
         test_  = "Kronecker product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= kron( eval( lhs_ ), eval( rhs_ ) );
            sres_   -= kron( eval( lhs_ ), eval( rhs_ ) );
            refres_ -= kron( eval( reflhs_ ), eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= kron( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   -= kron( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ -= kron( eval( treflhs_ ), eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Kronecker product with multiplication assignment
      //=====================================================================================

      // Kronecker product with multiplication assignment with the given vectors
      {
         test_  = "Kronecker product with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= kron( lhs_, rhs_ );
            sres_   *= kron( lhs_, rhs_ );
            refres_ *= kron( reflhs_, refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= kron( tlhs_, trhs_ );
            tsres_   *= kron( tlhs_, trhs_ );
            trefres_ *= kron( treflhs_, trefrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Kronecker product with multiplication assignment with evaluated vectors
      {
         test_  = "Kronecker product with multiplication assignment with evaluated vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= kron( eval( lhs_ ), eval( rhs_ ) );
            sres_   *= kron( eval( lhs_ ), eval( rhs_ ) );
            refres_ *= kron( eval( reflhs_ ), eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= kron( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   *= kron( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ *= kron( eval( treflhs_ ), eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the negated sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the negated vector Kronecker product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the Kronecker product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      //=====================================================================================
      // Negated Kronecker product
      //=====================================================================================

      // Negated Kronecker product with the given vectors
      {
         test_  = "Negated Kronecker product with the givven types";
         error_ = "Failed Kronecker product operation";

         try {
            initResults();
            dres_   = -kron( lhs_, rhs_ );
            sres_   = -kron( lhs_, rhs_ );
            refres_ = -kron( reflhs_, refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = -kron( tlhs_, trhs_ );
            tsres_   = -kron( tlhs_, trhs_ );
            trefres_ = -kron( treflhs_, trefrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Negated Kronecker product with evaluated vectors
      {
         test_  = "Negated Kronecker product with evaluated vectors";
         error_ = "Failed Kronecker product operation";

         try {
            initResults();
            dres_   = -kron( eval( lhs_ ), eval( rhs_ ) );
            sres_   = -kron( eval( lhs_ ), eval( rhs_ ) );
            refres_ = -kron( eval( reflhs_ ), eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = -kron( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   = -kron( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ = -kron( eval( treflhs_ ), eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Negated Kronecker product with addition assignment
      //=====================================================================================

      // Negated Kronecker product with addition assignment with the given vectors
      {
         test_  = "Negated Kronecker product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -kron( lhs_, rhs_ );
            sres_   += -kron( lhs_, rhs_ );
            refres_ += -kron( reflhs_, refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += -kron( tlhs_, trhs_ );
            tsres_   += -kron( tlhs_, trhs_ );
            trefres_ += -kron( treflhs_, trefrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Negated Kronecker product with addition assignment with evaluated vectors
      {
         test_  = "Negated Kronecker product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -kron( eval( lhs_ ), eval( rhs_ ) );
            sres_   += -kron( eval( lhs_ ), eval( rhs_ ) );
            refres_ += -kron( eval( reflhs_ ), eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += -kron( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   += -kron( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ += -kron( eval( treflhs_ ), eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Negated Kronecker product with subtraction assignment
      //=====================================================================================

      // Negated Kronecker product with subtraction assignment with the given vectors
      {
         test_  = "Negated Kronecker product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -kron( lhs_, rhs_ );
            sres_   -= -kron( lhs_, rhs_ );
            refres_ -= -kron( reflhs_, refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();
         try {
            initTransposeResults();
            tdres_   -= -kron( tlhs_, trhs_ );
            tsres_   -= -kron( tlhs_, trhs_ );
            trefres_ -= -kron( treflhs_, trefrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Negated Kronecker product with subtraction assignment with evaluated vectors
      {
         test_  = "Negated Kronecker product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -kron( eval( lhs_ ), eval( rhs_ ) );
            sres_   -= -kron( eval( lhs_ ), eval( rhs_ ) );
            refres_ -= -kron( eval( reflhs_ ), eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= -kron( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   -= -kron( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ -= -kron( eval( treflhs_ ), eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Negated Kronecker product with multiplication assignment
      //=====================================================================================

      // Negated Kronecker product with multiplication assignment with the given vectors
      {
         test_  = "Negated Kronecker product with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= -kron( lhs_, rhs_ );
            sres_   *= -kron( lhs_, rhs_ );
            refres_ *= -kron( reflhs_, refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= -kron( tlhs_, trhs_ );
            tsres_   *= -kron( tlhs_, trhs_ );
            trefres_ *= -kron( treflhs_, trefrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Negated Kronecker product with multiplication assignment with evaluated vectors
      {
         test_  = "Negated Kronecker product with multiplication assignment with evaluated vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= -kron( eval( lhs_ ), eval( rhs_ ) );
            sres_   *= -kron( eval( lhs_ ), eval( rhs_ ) );
            refres_ *= -kron( eval( reflhs_ ), eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= -kron( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   *= -kron( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ *= -kron( eval( treflhs_ ), eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the scaled sparse vector/dense vector Kronecker product.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the scaled vector Kronecker product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the Kronecker product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
template< typename T >    // Type of the scalar
void OperationTest<VT1,VT2>::testScaledOperation( T scalar )
{
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE( T );

   if( scalar == T(0) )
      throw std::invalid_argument( "Invalid scalar parameter" );


#if BLAZETEST_MATHTEST_TEST_SCALED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SCALED_OPERATION > 1 )
   {
      //=====================================================================================
      // Self-scaling (v*=s)
      //=====================================================================================

      // Self-scaling (v*=s)
      {
         test_ = "Self-scaling (v*=s)";

         try {
            dres_   = kron( lhs_, rhs_ );
            sres_   = dres_;
            refres_ = dres_;

            dres_   *= scalar;
            sres_   *= scalar;
            refres_ *= scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed self-scaling operation\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Scalar = " << scalar << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();
      }


      //=====================================================================================
      // Self-scaling (v=v*s)
      //=====================================================================================

      // Self-scaling (v=v*s)
      {
         test_ = "Self-scaling (v=v*s)";

         try {
            dres_   = kron( lhs_, rhs_ );
            sres_   = dres_;
            refres_ = dres_;

            dres_   = dres_   * scalar;
            sres_   = sres_   * scalar;
            refres_ = refres_ * scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed self-scaling operation\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Scalar = " << scalar << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();
      }


      //=====================================================================================
      // Self-scaling (v=s*v)
      //=====================================================================================

      // Self-scaling (v=s*v)
      {
         test_ = "Self-scaling (v=s*v)";

         try {
            dres_   = kron( lhs_, rhs_ );
            sres_   = dres_;
            refres_ = dres_;

            dres_   = scalar * dres_;
            sres_   = scalar * sres_;
            refres_ = scalar * refres_;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed self-scaling operation\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Scalar = " << scalar << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();
      }


      //=====================================================================================
      // Self-scaling (v/=s)
      //=====================================================================================

      // Self-scaling (v/=s)
      {
         test_ = "Self-scaling (v/=s)";

         try {
            dres_   = kron( lhs_, rhs_ );
            sres_   = dres_;
            refres_ = dres_;

            dres_   /= scalar;
            sres_   /= scalar;
            refres_ /= scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed self-scaling operation\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Scalar = " << scalar << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();
      }


      //=====================================================================================
      // Self-scaling (v=v/s)
      //=====================================================================================

      // Self-scaling (v=v/s)
      {
         test_ = "Self-scaling (v=v/s)";

         try {
            dres_   = kron( lhs_, rhs_ );
            sres_   = dres_;
            refres_ = dres_;

            dres_   = dres_   / scalar;
            sres_   = sres_   / scalar;
            refres_ = refres_ / scalar;
         }
         catch( std::exception& ex ) {
            std::ostringstream oss;
            oss << " Test : " << test_ << "\n"
                << " Error: Failed self-scaling operation\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Scalar = " << scalar << "\n"
                << "   Error message: " << ex.what() << "\n";
            throw std::runtime_error( oss.str() );
         }

         checkResults<VT1,VT2>();
      }


      //=====================================================================================
      // Scaled addition (s*OP)
      //=====================================================================================

      // Scaled addition with the given vectors
      {
         test_  = "Scaled addition with the given vectors (s*OP)";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = scalar * kron( lhs_, rhs_ );
            sres_   = scalar * kron( lhs_, rhs_ );
            refres_ = scalar * kron( reflhs_, refrhs_ );
         }
         catch( std::exception& ex ) {
             convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = scalar * kron( tlhs_, trhs_ );
            tsres_   = scalar * kron( tlhs_, trhs_ );
            trefres_ = scalar * kron( treflhs_, trefrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled addition with evaluated vectors
      {
         test_ = "Scaled addition with evaluated vectors (s*OP)";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = scalar * kron( eval( lhs_ ), eval( rhs_ ) );
            sres_   = scalar * kron( eval( lhs_ ), eval( rhs_ ) );
            refres_ = scalar * kron( eval( reflhs_ ), eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = scalar * kron( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   = scalar * kron( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ = scalar * kron( eval( treflhs_ ), eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled addition (OP*s)
      //=====================================================================================

      // Scaled addition with the given vectors
      {
         test_  = "Scaled addition with the given vectors (OP*s)";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = kron( lhs_, rhs_ ) * scalar;
            sres_   = kron( lhs_, rhs_ ) * scalar;
            refres_ = kron( reflhs_, refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = kron( tlhs_, trhs_ ) * scalar;
            tsres_   = kron( tlhs_, trhs_ ) * scalar;
            trefres_ = kron( treflhs_, trefrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled addition with evaluated vectors
      {
         test_  = "Scaled addition with evaluated vectors (OP*s)";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = kron( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            sres_   = kron( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            refres_ = kron( eval( reflhs_ ), eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = kron( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            tsres_   = kron( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            trefres_ = kron( eval( treflhs_ ), eval( trefrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled addition (OP/s)
      //=====================================================================================

      // Scaled addition with the given vectors
      {
         test_  = "Scaled addition with the given vectors (OP/s)";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = kron( lhs_, rhs_ ) / scalar;
            sres_   = kron( lhs_, rhs_ ) / scalar;
            refres_ = kron( reflhs_, refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = kron( tlhs_, trhs_ ) / scalar;
            tsres_   = kron( tlhs_, trhs_ ) / scalar;
            trefres_ = kron( treflhs_, trefrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled addition with evaluated vectors
      {
         test_  = "Scaled addition with evaluated vectors (OP/s)";
         error_ = "Failed addition operation";

         try {
            initResults();
            dres_   = kron( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            sres_   = kron( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            refres_ = kron( eval( reflhs_ ), eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = kron( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            tsres_   = kron( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            trefres_ = kron( eval( treflhs_ ), eval( trefrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled addition with addition assignment (s*OP)
      //=====================================================================================

      // Scaled addition with addition assignment with the given vectors
      {
         test_  = "Scaled addition with addition assignment with the given vectors (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * kron( lhs_, rhs_ );
            sres_   += scalar * kron( lhs_, rhs_ );
            refres_ += scalar * kron( reflhs_, refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += scalar * kron( tlhs_, trhs_ );
            tsres_   += scalar * kron( tlhs_, trhs_ );
            trefres_ += scalar * kron( treflhs_, trefrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled addition with addition assignment with evaluated vectors
      {
         test_  = "Scaled addition with addition assignment with evaluated vectors (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * kron( eval( lhs_ ), eval( rhs_ ) );
            sres_   += scalar * kron( eval( lhs_ ), eval( rhs_ ) );
            refres_ += scalar * kron( eval( reflhs_ ), eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += scalar * kron( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   += scalar * kron( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ += scalar * kron( eval( treflhs_ ), eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled addition with addition assignment (OP*s)
      //=====================================================================================

      // Scaled addition with addition assignment with the given vectors
      {
         test_  = "Scaled addition with addition assignment with the given vectors (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += kron( lhs_, rhs_ ) * scalar;
            sres_   += kron( lhs_, rhs_ ) * scalar;
            refres_ += kron( reflhs_, refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += kron( tlhs_, trhs_ ) * scalar;
            tsres_   += kron( tlhs_, trhs_ ) * scalar;
            trefres_ += kron( treflhs_, trefrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled addition with addition assignment with evaluated vectors
      {
         test_  = "Scaled addition with addition assignment with evaluated vectors (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += kron( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            sres_   += kron( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            refres_ += kron( eval( reflhs_ ), eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += kron( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            tsres_   += kron( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            trefres_ += kron( eval( treflhs_ ), eval( trefrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled addition with addition assignment (OP/s)
      //=====================================================================================

      // Scaled addition with addition assignment with the given vectors
      {
         test_  = "Scaled addition with addition assignment with the given vectors (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += kron( lhs_, rhs_ ) / scalar;
            sres_   += kron( lhs_, rhs_ ) / scalar;
            refres_ += kron( reflhs_, refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += kron( tlhs_, trhs_ ) / scalar;
            tsres_   += kron( tlhs_, trhs_ ) / scalar;
            trefres_ += kron( treflhs_, trefrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled addition with addition assignment with evaluated vectors
      {
         test_  = "Scaled addition with addition assignment with evaluated vectors (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += kron( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            sres_   += kron( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            refres_ += kron( eval( reflhs_ ), eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += kron( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            tsres_   += kron( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            trefres_ += kron( eval( treflhs_ ), eval( trefrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled addition with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled addition with subtraction assignment with the given vectors
      {
         test_  = "Scaled addition with subtraction assignment with the given vectors (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * kron( lhs_, rhs_ );
            sres_   -= scalar * kron( lhs_, rhs_ );
            refres_ -= scalar * kron( reflhs_, refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= scalar * kron( tlhs_, trhs_ );
            tsres_   -= scalar * kron( tlhs_, trhs_ );
            trefres_ -= scalar * kron( treflhs_, trefrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled addition with subtraction assignment with evaluated vectors
      {
         test_  = "Scaled addition with subtraction assignment with evaluated vectors (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * kron( eval( lhs_ ), eval( rhs_ ) );
            sres_   -= scalar * kron( eval( lhs_ ), eval( rhs_ ) );
            refres_ -= scalar * kron( eval( reflhs_ ), eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= scalar * kron( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   -= scalar * kron( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ -= scalar * kron( eval( treflhs_ ), eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled addition with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled addition with subtraction assignment with the given vectors
      {
         test_  = "Scaled addition with subtraction assignment with the given vectors (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= kron( lhs_, rhs_ ) * scalar;
            sres_   -= kron( lhs_, rhs_ ) * scalar;
            refres_ -= kron( reflhs_, refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= kron( tlhs_, trhs_ ) * scalar;
            tsres_   -= kron( tlhs_, trhs_ ) * scalar;
            trefres_ -= kron( treflhs_, trefrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled addition with subtraction assignment with evaluated vectors
      {
         test_  = "Scaled addition with subtraction assignment with evaluated vectors (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= kron( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            sres_   -= kron( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            refres_ -= kron( eval( reflhs_ ), eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= kron( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            tsres_   -= kron( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            trefres_ -= kron( eval( treflhs_ ), eval( trefrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled addition with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled addition with subtraction assignment with the given vectors
      {
         test_  = "Scaled addition with subtraction assignment with the given vectors (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= kron( lhs_, rhs_ ) / scalar;
            sres_   -= kron( lhs_, rhs_ ) / scalar;
            refres_ -= kron( reflhs_, refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= kron( tlhs_, trhs_ ) / scalar;
            tsres_   -= kron( tlhs_, trhs_ ) / scalar;
            trefres_ -= kron( treflhs_, trefrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled addition with subtraction assignment with evaluated vectors
      {
         test_  = "Scaled addition with subtraction assignment with evaluated vectors (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= kron( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            sres_   -= kron( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            refres_ -= kron( eval( reflhs_ ), eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= kron( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            tsres_   -= kron( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            trefres_ -= kron( eval( treflhs_ ), eval( trefrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled addition with multiplication assignment (s*OP)
      //=====================================================================================

      // Scaled addition with multiplication assignment with the given vectors
      {
         test_  = "Scaled addition with multiplication assignment with the given vectors (s*OP)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= scalar * kron( lhs_, rhs_ );
            sres_   *= scalar * kron( lhs_, rhs_ );
            refres_ *= scalar * kron( reflhs_, refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= scalar * kron( tlhs_, trhs_ );
            tsres_   *= scalar * kron( tlhs_, trhs_ );
            trefres_ *= scalar * kron( treflhs_, trefrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled addition with multiplication assignment with evaluated vectors
      {
         test_  = "Scaled addition with multiplication assignment with evaluated vectors (s*OP)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= scalar * kron( eval( lhs_ ), eval( rhs_ ) );
            sres_   *= scalar * kron( eval( lhs_ ), eval( rhs_ ) );
            refres_ *= scalar * kron( eval( reflhs_ ), eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= scalar * kron( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   *= scalar * kron( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ *= scalar * kron( eval( treflhs_ ), eval( trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled addition with multiplication assignment (OP*s)
      //=====================================================================================

      // Scaled addition with multiplication assignment with the given vectors
      {
         test_  = "Scaled addition with multiplication assignment with the given vectors (OP*s)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= kron( lhs_, rhs_ ) * scalar;
            sres_   *= kron( lhs_, rhs_ ) * scalar;
            refres_ *= kron( reflhs_, refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= kron( tlhs_, trhs_ ) * scalar;
            tsres_   *= kron( tlhs_, trhs_ ) * scalar;
            trefres_ *= kron( treflhs_, trefrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled addition with multiplication assignment with evaluated vectors
      {
         test_  = "Scaled addition with multiplication assignment with evaluated vectors (OP*s)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= kron( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            sres_   *= kron( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            refres_ *= kron( eval( reflhs_ ), eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= kron( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            tsres_   *= kron( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            trefres_ *= kron( eval( treflhs_ ), eval( trefrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled addition with multiplication assignment (OP/s)
      //=====================================================================================

      // Scaled addition with multiplication assignment with the given vectors
      {
         test_  = "Scaled addition with multiplication assignment with the given vectors (OP/s)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= kron( lhs_, rhs_ ) / scalar;
            sres_   *= kron( lhs_, rhs_ ) / scalar;
            refres_ *= kron( reflhs_, refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= kron( tlhs_, trhs_ ) / scalar;
            tsres_   *= kron( tlhs_, trhs_ ) / scalar;
            trefres_ *= kron( treflhs_, trefrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled addition with multiplication assignment with evaluated vectors
      {
         test_  = "Scaled addition with multiplication assignment with evaluated vectors (OP/s)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= kron( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            sres_   *= kron( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            refres_ *= kron( eval( reflhs_ ), eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= kron( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            tsres_   *= kron( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            trefres_ *= kron( eval( treflhs_ ), eval( trefrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the transpose sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the transpose vector Kronecker product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the Kronecker product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANS_OPERATION > 1 )
   {
      //=====================================================================================
      // Transpose Kronecker product
      //=====================================================================================

      // Transpose Kronecker product with the given vectors
      {
         test_  = "Transpose Kronecker product with the given vectors";
         error_ = "Failed Kronecker product operation";

         try {
            initTransposeResults();
            tdres_   = trans( kron( lhs_, rhs_ ) );
            tsres_   = trans( kron( lhs_, rhs_ ) );
            trefres_ = trans( kron( reflhs_, refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   = trans( kron( tlhs_, trhs_ ) );
            sres_   = trans( kron( tlhs_, trhs_ ) );
            refres_ = trans( kron( treflhs_, trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Transpose Kronecker product with evaluated vectors
      {
         test_  = "Transpose Kronecker product with evaluated vectors";
         error_ = "Failed Kronecker product operation";

         try {
            initTransposeResults();
            tdres_   = trans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   = trans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ = trans( kron( eval( reflhs_ ), eval( refrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   = trans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   = trans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ = trans( kron( eval( treflhs_ ), eval( trefrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Transpose Kronecker product with addition assignment
      //=====================================================================================

      // Transpose Kronecker product with addition assignment with the given vectors
      {
         test_  = "Transpose Kronecker product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initTransposeResults();
            tdres_   += trans( kron( lhs_, rhs_ ) );
            tsres_   += trans( kron( lhs_, rhs_ ) );
            trefres_ += trans( kron( reflhs_, refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   += trans( kron( tlhs_, trhs_ ) );
            sres_   += trans( kron( tlhs_, trhs_ ) );
            refres_ += trans( kron( treflhs_, trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Transpose Kronecker product with addition assignment with evaluated vectors
      {
         test_  = "Transpose Kronecker product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initTransposeResults();
            tdres_   += trans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   += trans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ += trans( kron( eval( reflhs_ ), eval( refrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   += trans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   += trans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ += trans( kron( eval( treflhs_ ), eval( trefrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Transpose Kronecker product with subtraction assignment
      //=====================================================================================

      // Transpose Kronecker product with subtraction assignment with the given vectors
      {
         test_  = "Transpose Kronecker product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initTransposeResults();
            tdres_   -= trans( kron( lhs_, rhs_ ) );
            tsres_   -= trans( kron( lhs_, rhs_ ) );
            trefres_ -= trans( kron( reflhs_, refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   -= trans( kron( tlhs_, trhs_ ) );
            sres_   -= trans( kron( tlhs_, trhs_ ) );
            refres_ -= trans( kron( treflhs_, trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Transpose Kronecker product with subtraction assignment with evaluated vectors
      {
         test_  = "Transpose Kronecker product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initTransposeResults();
            tdres_   -= trans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   -= trans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ -= trans( kron( eval( reflhs_ ), eval( refrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   -= trans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   -= trans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ -= trans( kron( eval( treflhs_ ), eval( trefrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Transpose Kronecker product with multiplication assignment
      //=====================================================================================

      // Transpose Kronecker product with multiplication assignment with the given vectors
      {
         test_  = "Transpose Kronecker product with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initTransposeResults();
            tdres_   *= trans( kron( lhs_, rhs_ ) );
            tsres_   *= trans( kron( lhs_, rhs_ ) );
            trefres_ *= trans( kron( reflhs_, refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   *= trans( kron( tlhs_, trhs_ ) );
            sres_   *= trans( kron( tlhs_, trhs_ ) );
            refres_ *= trans( kron( treflhs_, trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Transpose Kronecker product with multiplication assignment with evaluated vectors
      {
         test_  = "Transpose Kronecker product with multiplication assignment with evaluated vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initTransposeResults();
            tdres_   *= trans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   *= trans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ *= trans( kron( eval( reflhs_ ), eval( refrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   *= trans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   *= trans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ *= trans( kron( eval( treflhs_ ), eval( trefrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the conjugate transpose sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the conjugate transpose vector Kronecker product with plain assignment,
// addition assignment, subtraction assignment, and multiplication assignment. In case any
// error resulting from the Kronecker product or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testCTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION > 1 )
   {
      //=====================================================================================
      // Conjugate transpose Kronecker product
      //=====================================================================================

      // Conjugate transpose Kronecker product with the given vectors
      {
         test_  = "Conjugate transpose Kronecker product with the given vectors";
         error_ = "Failed Kronecker product operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( kron( lhs_, rhs_ ) );
            tsres_   = ctrans( kron( lhs_, rhs_ ) );
            trefres_ = ctrans( kron( reflhs_, refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   = ctrans( kron( tlhs_, trhs_ ) );
            sres_   = ctrans( kron( tlhs_, trhs_ ) );
            refres_ = ctrans( kron( treflhs_, trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Conjugate transpose Kronecker product with evaluated vectors
      {
         test_  = "Conjugate transpose Kronecker product with evaluated vectors";
         error_ = "Failed Kronecker product operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   = ctrans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ = ctrans( kron( eval( reflhs_ ), eval( refrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   = ctrans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   = ctrans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ = ctrans( kron( eval( treflhs_ ), eval( trefrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Conjugate transpose Kronecker product with addition assignment
      //=====================================================================================

      // Conjugate transpose Kronecker product with addition assignment with the given vectors
      {
         test_  = "Conjugate transpose Kronecker product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initTransposeResults();
            tdres_   += ctrans( kron( lhs_, rhs_ ) );
            tsres_   += ctrans( kron( lhs_, rhs_ ) );
            trefres_ += ctrans( kron( reflhs_, refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   += ctrans( kron( tlhs_, trhs_ ) );
            sres_   += ctrans( kron( tlhs_, trhs_ ) );
            refres_ += ctrans( kron( treflhs_, trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Conjugate transpose Kronecker product with addition assignment with evaluated vectors
      {
         test_  = "Conjugate transpose Kronecker product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initTransposeResults();
            tdres_   += ctrans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   += ctrans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ += ctrans( kron( eval( reflhs_ ), eval( refrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   += ctrans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   += ctrans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ += ctrans( kron( eval( treflhs_ ), eval( trefrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Conjugate transpose Kronecker product with subtraction assignment
      //=====================================================================================

      // Conjugate transpose Kronecker product with subtraction assignment with the given vectors
      {
         test_  = "Conjugate transpose Kronecker product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initTransposeResults();
            tdres_   -= ctrans( kron( lhs_, rhs_ ) );
            tsres_   -= ctrans( kron( lhs_, rhs_ ) );
            trefres_ -= ctrans( kron( reflhs_, refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   -= ctrans( kron( tlhs_, trhs_ ) );
            sres_   -= ctrans( kron( tlhs_, trhs_ ) );
            refres_ -= ctrans( kron( treflhs_, trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Conjugate transpose Kronecker product with subtraction assignment with evaluated vectors
      {
         test_  = "Conjugate transpose Kronecker product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initTransposeResults();
            tdres_   -= ctrans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   -= ctrans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ -= ctrans( kron( eval( reflhs_ ), eval( refrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   -= ctrans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   -= ctrans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ -= ctrans( kron( eval( treflhs_ ), eval( trefrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Conjugate transpose Kronecker product with multiplication assignment
      //=====================================================================================

      // Conjugate transpose Kronecker product with multiplication assignment with the given vectors
      {
         test_  = "Conjugate transpose Kronecker product with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initTransposeResults();
            tdres_   *= ctrans( kron( lhs_, rhs_ ) );
            tsres_   *= ctrans( kron( lhs_, rhs_ ) );
            trefres_ *= ctrans( kron( reflhs_, refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   *= ctrans( kron( tlhs_, trhs_ ) );
            sres_   *= ctrans( kron( tlhs_, trhs_ ) );
            refres_ *= ctrans( kron( treflhs_, trefrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Conjugate transpose Kronecker product with multiplication assignment with evaluated vectors
      {
         test_  = "Conjugate transpose Kronecker product with multiplication assignment with evaluated vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initTransposeResults();
            tdres_   *= ctrans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   *= ctrans( kron( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ *= ctrans( kron( eval( reflhs_ ), eval( refrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   *= ctrans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   *= ctrans( kron( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ *= ctrans( kron( eval( treflhs_ ), eval( trefrhs_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the abs sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the abs vector Kronecker product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the Kronecker product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testAbsOperation()
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
/*!\brief Testing the conjugate sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the conjugate vector Kronecker product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the Kronecker product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testConjOperation()
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
/*!\brief Testing the \a real sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the \a real vector Kronecker product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the Kronecker product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testRealOperation()
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
/*!\brief Testing the \a imag sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the \a imag vector Kronecker product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the Kronecker product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testImagOperation()
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
/*!\brief Testing the evaluated sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the evaluated vector Kronecker product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the Kronecker product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testEvalOperation()
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
/*!\brief Testing the serialized sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the serialized vector Kronecker product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the Kronecker product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testSerialOperation()
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
/*!\brief Testing the non-aliased sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the non-aliased vector Kronecker product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the Kronecker product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testNoAliasOperation()
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
/*!\brief Testing the non-SIMD sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the non-SIMD vector Kronecker product with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the Kronecker product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testNoSIMDOperation()
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
/*!\brief Testing the subvector-wise sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the subvector-wise vector Kronecker product with plain assignment,
// addition assignment, subtraction assignment, and multiplication assignment. In case any
// error resulting from the Kronecker product or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testSubvectorOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION > 1 )
   {
      if( lhs_.size() * rhs_.size() == 0UL )
         return;


      //=====================================================================================
      // Subvector-wise Kronecker product
      //=====================================================================================

      // Subvector-wise Kronecker product with the given vectors
      {
         test_  = "Subvector-wise Kronecker product with the given vectors";
         error_ = "Failed Kronecker product operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size()*rhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size()*rhs_.size() - index );
               subvector( dres_  , index, size ) = subvector( kron( lhs_, rhs_ )      , index, size );
               subvector( sres_  , index, size ) = subvector( kron( lhs_, rhs_ )      , index, size );
               subvector( refres_, index, size ) = subvector( kron( reflhs_, refrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size()*trhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size()*trhs_.size() - index );
               subvector( tdres_  , index, size ) = subvector( kron( tlhs_, trhs_ )      , index, size );
               subvector( tsres_  , index, size ) = subvector( kron( tlhs_, trhs_ )      , index, size );
               subvector( trefres_, index, size ) = subvector( kron( treflhs_, trefrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subvector-wise Kronecker product with evaluated vectors
      {
         test_  = "Subvector-wise Kronecker product with evaluated vectors";
         error_ = "Failed Kronecker product operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size()*rhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size()*rhs_.size() - index );
               subvector( dres_  , index, size ) = subvector( kron( eval( lhs_ ), eval( rhs_ ) )      , index, size );
               subvector( sres_  , index, size ) = subvector( kron( eval( lhs_ ), eval( rhs_ ) )      , index, size );
               subvector( refres_, index, size ) = subvector( kron( eval( reflhs_ ), eval( refrhs_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size()*trhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size()*trhs_.size() - index );
               subvector( tdres_  , index, size ) = subvector( kron( eval( tlhs_ ), eval( trhs_ ) )      , index, size );
               subvector( tsres_  , index, size ) = subvector( kron( eval( tlhs_ ), eval( trhs_ ) )      , index, size );
               subvector( trefres_, index, size ) = subvector( kron( eval( treflhs_ ), eval( trefrhs_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Subvector-wise Kronecker product with addition assignment
      //=====================================================================================

      // Subvector-wise Kronecker product with addition assignment with the given vectors
      {
         test_  = "Subvector-wise Kronecker product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size()*rhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size()*rhs_.size() - index );
               subvector( dres_  , index, size ) += subvector( kron( lhs_, rhs_ )      , index, size );
               subvector( sres_  , index, size ) += subvector( kron( lhs_, rhs_ )      , index, size );
               subvector( refres_, index, size ) += subvector( kron( reflhs_, refrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size()*trhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size()*trhs_.size() - index );
               subvector( tdres_  , index, size ) += subvector( kron( tlhs_, trhs_ )      , index, size );
               subvector( tsres_  , index, size ) += subvector( kron( tlhs_, trhs_ )      , index, size );
               subvector( trefres_, index, size ) += subvector( kron( treflhs_, trefrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subvector-wise Kronecker product with addition assignment with evaluated vectors
      {
         test_  = "Subvector-wise Kronecker product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size()*rhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size()*rhs_.size() - index );
               subvector( dres_  , index, size ) += subvector( kron( eval( lhs_ ), eval( rhs_ ) )      , index, size );
               subvector( sres_  , index, size ) += subvector( kron( eval( lhs_ ), eval( rhs_ ) )      , index, size );
               subvector( refres_, index, size ) += subvector( kron( eval( reflhs_ ), eval( refrhs_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size()*trhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size()*trhs_.size() - index );
               subvector( tdres_  , index, size ) += subvector( kron( eval( tlhs_ ), eval( trhs_ ) )      , index, size );
               subvector( tsres_  , index, size ) += subvector( kron( eval( tlhs_ ), eval( trhs_ ) )      , index, size );
               subvector( trefres_, index, size ) += subvector( kron( eval( treflhs_ ), eval( trefrhs_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Subvector-wise Kronecker product with subtraction assignment
      //=====================================================================================

      // Subvector-wise Kronecker product with subtraction assignment with the given vectors
      {
         test_  = "Subvector-wise Kronecker product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size()*rhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size()*rhs_.size() - index );
               subvector( dres_  , index, size ) -= subvector( kron( lhs_, rhs_ )      , index, size );
               subvector( sres_  , index, size ) -= subvector( kron( lhs_, rhs_ )      , index, size );
               subvector( refres_, index, size ) -= subvector( kron( reflhs_, refrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size()*trhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size()*trhs_.size() - index );
               subvector( tdres_  , index, size ) -= subvector( kron( tlhs_, trhs_ )      , index, size );
               subvector( tsres_  , index, size ) -= subvector( kron( tlhs_, trhs_ )      , index, size );
               subvector( trefres_, index, size ) -= subvector( kron( treflhs_, trefrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subvector-wise Kronecker product with subtraction assignment with evaluated vectors
      {
         test_  = "Subvector-wise Kronecker product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size()*rhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size()*rhs_.size() - index );
               subvector( dres_  , index, size ) -= subvector( kron( eval( lhs_ ), eval( rhs_ ) )      , index, size );
               subvector( sres_  , index, size ) -= subvector( kron( eval( lhs_ ), eval( rhs_ ) )      , index, size );
               subvector( refres_, index, size ) -= subvector( kron( eval( reflhs_ ), eval( refrhs_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size()*trhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size()*trhs_.size() - index );
               subvector( tdres_  , index, size ) -= subvector( kron( eval( tlhs_ ), eval( trhs_ ) )      , index, size );
               subvector( tsres_  , index, size ) -= subvector( kron( eval( tlhs_ ), eval( trhs_ ) )      , index, size );
               subvector( trefres_, index, size ) -= subvector( kron( eval( treflhs_ ), eval( trefrhs_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Subvector-wise Kronecker product with multiplication assignment
      //=====================================================================================

      // Subvector-wise Kronecker product with multiplication assignment with the given vectors
      {
         test_  = "Subvector-wise Kronecker product with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size()*rhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size()*rhs_.size() - index );
               subvector( dres_  , index, size ) *= subvector( kron( lhs_, rhs_ )      , index, size );
               subvector( sres_  , index, size ) *= subvector( kron( lhs_, rhs_ )      , index, size );
               subvector( refres_, index, size ) *= subvector( kron( reflhs_, refrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size()*trhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size()*trhs_.size() - index );
               subvector( tdres_  , index, size ) *= subvector( kron( tlhs_, trhs_ )      , index, size );
               subvector( tsres_  , index, size ) *= subvector( kron( tlhs_, trhs_ )      , index, size );
               subvector( trefres_, index, size ) *= subvector( kron( treflhs_, trefrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subvector-wise Kronecker product with multiplication assignment with evaluated vectors
      {
         test_  = "Subvector-wise Kronecker product with multiplication assignment with evaluated vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size()*rhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size()*rhs_.size() - index );
               subvector( dres_  , index, size ) *= subvector( kron( eval( lhs_ ), eval( rhs_ ) )      , index, size );
               subvector( sres_  , index, size ) *= subvector( kron( eval( lhs_ ), eval( rhs_ ) )      , index, size );
               subvector( refres_, index, size ) *= subvector( kron( eval( reflhs_ ), eval( refrhs_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size()*trhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size()*trhs_.size() - index );
               subvector( tdres_  , index, size ) *= subvector( kron( eval( tlhs_ ), eval( trhs_ ) )      , index, size );
               subvector( tsres_  , index, size ) *= subvector( kron( eval( tlhs_ ), eval( trhs_ ) )      , index, size );
               subvector( trefres_, index, size ) *= subvector( kron( eval( treflhs_ ), eval( trefrhs_ ) ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      try {
         auto sv = subvector( kron( lhs_, rhs_ ), 1UL, lhs_.size() * rhs_.size() );

         std::ostringstream oss;
         oss << " Test: Subvector construction\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n"
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
         auto sv = subvector( kron( lhs_, rhs_ ), lhs_.size() * rhs_.size(), 1UL );

         std::ostringstream oss;
         oss << " Test: Subvector construction\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n"
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
/*!\brief Skipping the subvector-wise sparse vector/dense vector Kronecker product.
//
// \return void
//
// This function is called in case the subvector-wise vector/vector Kronecker product operation
// is not available for the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testSubvectorOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the elements-wise sparse vector/dense vector Kronecker product.
//
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the elements-wise vector Kronecker product with plain assignment,
// addition assignment, subtraction assignment, and multiplication assignment. In case any
// error resulting from the Kronecker product or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testElementsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_ELEMENTS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ELEMENTS_OPERATION > 1 )
   {
      if( lhs_.size() == 0UL || rhs_.size() == 0UL )
         return;


      std::vector<size_t> indices( lhs_.size() * rhs_.size() );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );


      //=====================================================================================
      // Elements-wise Kronecker product
      //=====================================================================================

      // Elements-wise Kronecker product with the given vectors
      {
         test_  = "Elements-wise Kronecker product with the given vectors";
         error_ = "Failed Kronecker product operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) = elements( kron( lhs_, rhs_ )      , &indices[index], n );
               elements( sres_  , &indices[index], n ) = elements( kron( lhs_, rhs_ )      , &indices[index], n );
               elements( refres_, &indices[index], n ) = elements( kron( reflhs_, refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) = elements( kron( tlhs_, trhs_ )      , &indices[index], n );
               elements( tsres_  , &indices[index], n ) = elements( kron( tlhs_, trhs_ )      , &indices[index], n );
               elements( trefres_, &indices[index], n ) = elements( kron( treflhs_, trefrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Elements-wise Kronecker product with evaluated vectors
      {
         test_  = "Elements-wise Kronecker product with evaluated vectors";
         error_ = "Failed Kronecker product operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) = elements( kron( eval( lhs_ ), eval( rhs_ ) )      , &indices[index], n );
               elements( sres_  , &indices[index], n ) = elements( kron( eval( lhs_ ), eval( rhs_ ) )      , &indices[index], n );
               elements( refres_, &indices[index], n ) = elements( kron( eval( reflhs_ ), eval( refrhs_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) = elements( kron( eval( tlhs_ ), eval( trhs_ ) )      , &indices[index], n );
               elements( tsres_  , &indices[index], n ) = elements( kron( eval( tlhs_ ), eval( trhs_ ) )      , &indices[index], n );
               elements( trefres_, &indices[index], n ) = elements( kron( eval( treflhs_ ), eval( trefrhs_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Elements-wise Kronecker product with addition assignment
      //=====================================================================================

      // Elements-wise Kronecker product with addition assignment with the given vectors
      {
         test_  = "Elements-wise Kronecker product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) += elements( kron( lhs_, rhs_ )      , &indices[index], n );
               elements( sres_  , &indices[index], n ) += elements( kron( lhs_, rhs_ )      , &indices[index], n );
               elements( refres_, &indices[index], n ) += elements( kron( reflhs_, refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) += elements( kron( tlhs_, trhs_ )      , &indices[index], n );
               elements( tsres_  , &indices[index], n ) += elements( kron( tlhs_, trhs_ )      , &indices[index], n );
               elements( trefres_, &indices[index], n ) += elements( kron( treflhs_, trefrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Elements-wise Kronecker product with addition assignment with evaluated vectors
      {
         test_  = "Elements-wise Kronecker product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) += elements( kron( eval( lhs_ ), eval( rhs_ ) )      , &indices[index], n );
               elements( sres_  , &indices[index], n ) += elements( kron( eval( lhs_ ), eval( rhs_ ) )      , &indices[index], n );
               elements( refres_, &indices[index], n ) += elements( kron( eval( reflhs_ ), eval( refrhs_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) += elements( kron( eval( tlhs_ ), eval( trhs_ ) )      , &indices[index], n );
               elements( tsres_  , &indices[index], n ) += elements( kron( eval( tlhs_ ), eval( trhs_ ) )      , &indices[index], n );
               elements( trefres_, &indices[index], n ) += elements( kron( eval( treflhs_ ), eval( trefrhs_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Elements-wise Kronecker product with subtraction assignment
      //=====================================================================================

      // Elements-wise Kronecker product with subtraction assignment with the given vectors
      {
         test_  = "Elements-wise Kronecker product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) -= elements( kron( lhs_, rhs_ )      , &indices[index], n );
               elements( sres_  , &indices[index], n ) -= elements( kron( lhs_, rhs_ )      , &indices[index], n );
               elements( refres_, &indices[index], n ) -= elements( kron( reflhs_, refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) -= elements( kron( tlhs_, trhs_ )      , &indices[index], n );
               elements( tsres_  , &indices[index], n ) -= elements( kron( tlhs_, trhs_ )      , &indices[index], n );
               elements( trefres_, &indices[index], n ) -= elements( kron( treflhs_, trefrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Elements-wise Kronecker product with subtraction assignment with evaluated vectors
      {
         test_  = "Elements-wise Kronecker product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) -= elements( kron( eval( lhs_ ), eval( rhs_ ) )      , &indices[index], n );
               elements( sres_  , &indices[index], n ) -= elements( kron( eval( lhs_ ), eval( rhs_ ) )      , &indices[index], n );
               elements( refres_, &indices[index], n ) -= elements( kron( eval( reflhs_ ), eval( refrhs_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) -= elements( kron( eval( tlhs_ ), eval( trhs_ ) )      , &indices[index], n );
               elements( tsres_  , &indices[index], n ) -= elements( kron( eval( tlhs_ ), eval( trhs_ ) )      , &indices[index], n );
               elements( trefres_, &indices[index], n ) -= elements( kron( eval( treflhs_ ), eval( trefrhs_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Elements-wise Kronecker product with multiplication assignment
      //=====================================================================================

      // Elements-wise Kronecker product with multiplication assignment with the given vectors
      {
         test_  = "Elements-wise Kronecker product with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) *= elements( kron( lhs_, rhs_ )      , &indices[index], n );
               elements( sres_  , &indices[index], n ) *= elements( kron( lhs_, rhs_ )      , &indices[index], n );
               elements( refres_, &indices[index], n ) *= elements( kron( reflhs_, refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) *= elements( kron( tlhs_, trhs_ )      , &indices[index], n );
               elements( tsres_  , &indices[index], n ) *= elements( kron( tlhs_, trhs_ )      , &indices[index], n );
               elements( trefres_, &indices[index], n ) *= elements( kron( treflhs_, trefrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Elements-wise Kronecker product with multiplication assignment with evaluated vectors
      {
         test_  = "Elements-wise Kronecker product with multiplication assignment with evaluated vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) *= elements( kron( eval( lhs_ ), eval( rhs_ ) )      , &indices[index], n );
               elements( sres_  , &indices[index], n ) *= elements( kron( eval( lhs_ ), eval( rhs_ ) )      , &indices[index], n );
               elements( refres_, &indices[index], n ) *= elements( kron( eval( reflhs_ ), eval( refrhs_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( tdres_  , &indices[index], n ) *= elements( kron( eval( tlhs_ ), eval( trhs_ ) )      , &indices[index], n );
               elements( tsres_  , &indices[index], n ) *= elements( kron( eval( tlhs_ ), eval( trhs_ ) )      , &indices[index], n );
               elements( trefres_, &indices[index], n ) *= elements( kron( eval( treflhs_ ), eval( trefrhs_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      try {
         auto e = elements( kron( lhs_, rhs_ ), blaze::index_sequence<31*32>() );

         std::ostringstream oss;
         oss << " Test: Elements construction (index_sequence)\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n"
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
         auto e = elements( kron( lhs_, rhs_ ), { lhs_.size() * rhs_.size() } );

         std::ostringstream oss;
         oss << " Test: Elements construction (initializer_list)\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n"
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
         auto e = elements( kron( lhs_, rhs_ ), [index=lhs_.size()*rhs_.size()]( size_t ){ return index; }, 1UL );

         std::ostringstream oss;
         oss << " Test: Elements construction (lambda)\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side sparse vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n"
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
/*!\brief Skipping the elements-wise sparse vector/dense vector Kronecker product.
//
// \return void
//
// This function is called in case the elements-wise vector/vector Kronecker product operation
// is not available for the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testElementsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the customized sparse vector/dense vector Kronecker product.
//
// \param op The custom operation to be tested.
// \param name The human-readable name of the operation.
// \return void
// \exception std::runtime_error Kronecker product error detected.
//
// This function tests the vector Kronecker product with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment in combination with a custom operation.
// In case any error resulting from the Kronecker product or the subsequent assignment is
// detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
template< typename OP >   // Type of the custom operation
void OperationTest<VT1,VT2>::testCustomOperation( OP op, const std::string& name )
{
   //=====================================================================================
   // Customized Kronecker product
   //=====================================================================================

   // Customized Kronecker product with the given vectors
   {
      test_  = "Customized Kronecker product with the given vectors (" + name + ")";
      error_ = "Failed addition operation";

      try {
         initResults();
         dres_   = op( kron( lhs_, rhs_ ) );
         sres_   = op( kron( lhs_, rhs_ ) );
         refres_ = op( kron( reflhs_, refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   = op( kron( tlhs_, trhs_ ) );
         tsres_   = op( kron( tlhs_, trhs_ ) );
         trefres_ = op( kron( treflhs_, trefrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }

   // Customized Kronecker product with evaluated vectors
   {
      test_  = "Customized Kronecker product with evaluated vectors (" + name + ")";
      error_ = "Failed addition operation";

      try {
         initResults();
         dres_   = op( kron( eval( lhs_ ), eval( rhs_ ) ) );
         sres_   = op( kron( eval( lhs_ ), eval( rhs_ ) ) );
         refres_ = op( kron( eval( reflhs_ ), eval( refrhs_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   = op( kron( eval( tlhs_ ), eval( trhs_ ) ) );
         tsres_   = op( kron( eval( tlhs_ ), eval( trhs_ ) ) );
         trefres_ = op( kron( eval( treflhs_ ), eval( trefrhs_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }


   //=====================================================================================
   // Customized Kronecker product with addition assignment
   //=====================================================================================

   // Customized Kronecker product with addition assignment with the given vectors
   {
      test_  = "Customized Kronecker product with addition assignment with the given vectors (" + name + ")";
      error_ = "Failed addition assignment operation";

      try {
         initResults();
         dres_   += op( kron( lhs_, rhs_ ) );
         sres_   += op( kron( lhs_, rhs_ ) );
         refres_ += op( kron( reflhs_, refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   += op( kron( tlhs_, trhs_ ) );
         tsres_   += op( kron( tlhs_, trhs_ ) );
         trefres_ += op( kron( treflhs_, trefrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }

   // Customized Kronecker product with addition assignment with evaluated vectors
   {
      test_  = "Customized Kronecker product with addition assignment with evaluated vectors (" + name + ")";
      error_ = "Failed addition assignment operation";

      try {
         initResults();
         dres_   += op( kron( eval( lhs_ ), eval( rhs_ ) ) );
         sres_   += op( kron( eval( lhs_ ), eval( rhs_ ) ) );
         refres_ += op( kron( eval( reflhs_ ), eval( refrhs_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   += op( kron( eval( tlhs_ ), eval( trhs_ ) ) );
         tsres_   += op( kron( eval( tlhs_ ), eval( trhs_ ) ) );
         trefres_ += op( kron( eval( treflhs_ ), eval( trefrhs_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }


   //=====================================================================================
   // Customized Kronecker product with subtraction assignment
   //=====================================================================================

   // Customized Kronecker product with subtraction assignment with the given vectors
   {
      test_  = "Customized Kronecker product with subtraction assignment with the given vectors (" + name + ")";
      error_ = "Failed subtraction assignment operation";

      try {
         initResults();
         dres_   -= op( kron( lhs_, rhs_ ) );
         sres_   -= op( kron( lhs_, rhs_ ) );
         refres_ -= op( kron( reflhs_, refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   -= op( kron( tlhs_, trhs_ ) );
         tsres_   -= op( kron( tlhs_, trhs_ ) );
         trefres_ -= op( kron( treflhs_, trefrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }

   // Customized Kronecker product with subtraction assignment with evaluated vectors
   {
      test_  = "Customized Kronecker product with subtraction assignment with evaluated vectors (" + name + ")";
      error_ = "Failed subtraction assignment operation";

      try {
         initResults();
         dres_   -= op( kron( eval( lhs_ ), eval( rhs_ ) ) );
         sres_   -= op( kron( eval( lhs_ ), eval( rhs_ ) ) );
         refres_ -= op( kron( eval( reflhs_ ), eval( refrhs_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   -= op( kron( eval( tlhs_ ), eval( trhs_ ) ) );
         tsres_   -= op( kron( eval( tlhs_ ), eval( trhs_ ) ) );
         trefres_ -= op( kron( eval( treflhs_ ), eval( trefrhs_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }


   //=====================================================================================
   // Customized Kronecker product with multiplication assignment
   //=====================================================================================

   // Customized Kronecker product with multiplication assignment with the given vectors
   {
      test_  = "Customized Kronecker product with multiplication assignment with the given vectors (" + name + ")";
      error_ = "Failed multiplication assignment operation";

      try {
         initResults();
         dres_   *= op( kron( lhs_, rhs_ ) );
         sres_   *= op( kron( lhs_, rhs_ ) );
         refres_ *= op( kron( reflhs_, refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   *= op( kron( tlhs_, trhs_ ) );
         tsres_   *= op( kron( tlhs_, trhs_ ) );
         trefres_ *= op( kron( treflhs_, trefrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }

   // Customized Kronecker product with multiplication assignment with evaluated vectors
   {
      test_  = "Customized Kronecker product with multiplication assignment with evaluated vectors (" + name + ")";
      error_ = "Failed multiplication assignment operation";

      try {
         initResults();
         dres_   *= op( kron( eval( lhs_ ), eval( rhs_ ) ) );
         sres_   *= op( kron( eval( lhs_ ), eval( rhs_ ) ) );
         refres_ *= op( kron( eval( reflhs_ ), eval( refrhs_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   *= op( kron( eval( tlhs_ ), eval( trhs_ ) ) );
         tsres_   *= op( kron( eval( tlhs_ ), eval( trhs_ ) ) );
         trefres_ *= op( kron( eval( treflhs_ ), eval( trefrhs_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
template< typename LT     // Type of the left-hand side operand
        , typename RT >   // Type of the right-hand side operand
void OperationTest<VT1,VT2>::checkResults()
{
   using blaze::IsRowVector;

   if( !isEqual( dres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result vector detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side sparse " << ( IsRowVector<LT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side dense " << ( IsRowVector<RT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( RT ).name() << "\n"
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
          << "   Left-hand side sparse " << ( IsRowVector<LT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side dense " << ( IsRowVector<RT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( RT ).name() << "\n"
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
template< typename LT     // Type of the left-hand side operand
        , typename RT >   // Type of the right-hand side operand
void OperationTest<VT1,VT2>::checkTransposeResults()
{
   using blaze::IsRowVector;

   if( !isEqual( tdres_, trefres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result vector detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side sparse " << ( IsRowVector<LT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side dense " << ( IsRowVector<RT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Result:\n" << tdres_ << "\n"
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
          << "   Left-hand side sparse " << ( IsRowVector<LT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side dense " << ( IsRowVector<RT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Result:\n" << tsres_ << "\n"
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::initResults()
{
   const blaze::UnderlyingBuiltin_t<DRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<DRE> max( randmax );

   resize( dres_, size( lhs_ ) * size( rhs_ ) );
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::initTransposeResults()
{
   const blaze::UnderlyingBuiltin_t<TDRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<TDRE> max( randmax );

   resize( tdres_, size( tlhs_ ) * size( trhs_ ) );
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
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
template< typename LT     // Type of the left-hand side operand
        , typename RT >   // Type of the right-hand side operand
void OperationTest<VT1,VT2>::convertException( const std::exception& ex )
{
   using blaze::IsRowVector;

   std::ostringstream oss;
   oss << " Test : " << test_ << "\n"
       << " Error: " << error_ << "\n"
       << " Details:\n"
       << "   Random seed = " << blaze::getSeed() << "\n"
       << "   Left-hand side sparse " << ( IsRowVector<LT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
       << "     " << typeid( LT ).name() << "\n"
       << "   Right-hand side dense " << ( IsRowVector<RT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
       << "     " << typeid( RT ).name() << "\n"
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
/*!\brief Testing the vector Kronecker product between two specific vector types.
//
// \param creator1 The creator for the left-hand side sparse vector.
// \param creator2 The creator for the right-hand side dense vector.
// \return void
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , typename VT2 >  // Type of the right-hand side dense vector
void runTest( const Creator<VT1>& creator1, const Creator<VT2>& creator2 )
{
#if BLAZETEST_MATHTEST_TEST_MULTIPLICATION
   if( BLAZETEST_MATHTEST_TEST_MULTIPLICATION > 1 )
   {
      for( size_t rep=0UL; rep<BLAZETEST_REPETITIONS; ++rep ) {
         OperationTest<VT1,VT2>( creator1, creator2 );
      }
   }
#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the definition of a sparse vector/dense vector Kronecker product test case.
*/
#define DEFINE_SVECDVECKRON_OPERATION_TEST( VT1, VT2 ) \
   extern template class blazetest::mathtest::operations::svecdveckron::OperationTest<VT1,VT2>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a sparse vector/dense vector Kronecker product test case.
*/
#define RUN_SVECDVECKRON_OPERATION_TEST( C1, C2 ) \
   blazetest::mathtest::operations::svecdveckron::runTest( C1, C2 )
/*! \endcond */
//*************************************************************************************************

} // namespace svecdveckron

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
