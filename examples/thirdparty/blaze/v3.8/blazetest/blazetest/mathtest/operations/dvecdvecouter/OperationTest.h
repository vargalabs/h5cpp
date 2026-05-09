//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/dvecdvecouter/OperationTest.h
//  \brief Header file for the dense vector/dense vector outer product operation test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_DVECDVECOUTER_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_DVECDVECOUTER_OPERATIONTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>
#include <blaze/math/Aliases.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/Functors.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingScalar.h>
#include <blaze/math/Views.h>
#include <blaze/util/constraints/DerivedFrom.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/Not.h>
#include <blaze/util/mpl/Or.h>
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

namespace dvecdvecouter {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the dense vector/dense vector multiplication operation test.
//
// This class template represents one particular outer product test between two vectors of a
// particular type. The two template arguments \a VT1 and \a VT2 represent the types of the
// left-hand side and right-hand side vector, respectively.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
class OperationTest
{
 private:
   //**Type definitions****************************************************************************
   using ET1 = blaze::ElementType_t<VT1>;  //!< Element type 1
   using ET2 = blaze::ElementType_t<VT2>;  //!< Element type 2

   using TVT1 = blaze::TransposeType_t<VT1>;  //!< Transpose vector type 1
   using TVT2 = blaze::TransposeType_t<VT2>;  //!< Transpose vector type 2

   //! Dense result type
   using DRE = blaze::MultTrait_t<VT1,TVT2>;

   using DET   = blaze::ElementType_t<DRE>;     //!< Element type of the dense result
   using ODRE  = blaze::OppositeType_t<DRE>;    //!< Dense result type with opposite storage order
   using TDRE  = blaze::TransposeType_t<DRE>;   //!< Transpose dense result type
   using TODRE = blaze::TransposeType_t<ODRE>;  //!< Transpose dense result type with opposite storage order

   //! Sparse result type
   using SRE = blaze::CompressedMatrix<DET,false>;

   using SET   = blaze::ElementType_t<SRE>;     //!< Element type of the sparse result
   using OSRE  = blaze::OppositeType_t<SRE>;    //!< Sparse result type with opposite storage order
   using TSRE  = blaze::TransposeType_t<SRE>;   //!< Transpose sparse result type
   using TOSRE = blaze::TransposeType_t<OSRE>;  //!< Transpose sparse result type with opposite storage order

   using RT1 = blaze::DynamicVector<ET1,false>;    //!< Reference type 1
   using RT2 = blaze::CompressedVector<ET2,true>;  //!< Reference type 2
   using RRE = blaze::MultTrait_t<RT1,RT2>;        //!< Reference result type

   //! Type of the outer product expression
   using OuterExprType =
      blaze::RemoveCVRef_t< decltype( std::declval<VT1>() * std::declval<TVT2>() ) >;
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
                          void testDeclSymOperation  ( blaze::TrueType  );
                          void testDeclSymOperation  ( blaze::FalseType );
                          void testDeclHermOperation ( blaze::TrueType  );
                          void testDeclHermOperation ( blaze::FalseType );
                          void testDeclLowOperation  ( blaze::TrueType  );
                          void testDeclLowOperation  ( blaze::FalseType );
                          void testDeclUppOperation  ( blaze::TrueType  );
                          void testDeclUppOperation  ( blaze::FalseType );
                          void testDeclDiagOperation ( blaze::TrueType  );
                          void testDeclDiagOperation ( blaze::FalseType );
                          void testSubmatrixOperation( blaze::TrueType  );
                          void testSubmatrixOperation( blaze::FalseType );
                          void testRowOperation      ( blaze::TrueType  );
                          void testRowOperation      ( blaze::FalseType );
                          void testRowsOperation     ( blaze::TrueType  );
                          void testRowsOperation     ( blaze::FalseType );
                          void testColumnOperation   ( blaze::TrueType  );
                          void testColumnOperation   ( blaze::FalseType );
                          void testColumnsOperation  ( blaze::TrueType  );
                          void testColumnsOperation  ( blaze::FalseType );
                          void testBandOperation     ( blaze::TrueType  );
                          void testBandOperation     ( blaze::FalseType );

   template< typename OP > void testCustomOperation( OP op, const std::string& name );
   //@}
   //**********************************************************************************************

   //**Error detection functions*******************************************************************
   /*!\name Error detection functions */
   //@{
   void checkResults();
   void checkTransposeResults();
   void checkExceptionMessage( const std::exception& ex, const std::string& message );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initResults();
   void initTransposeResults();
   void convertException( const std::exception& ex );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VT1   lhs_;     //!< The left-hand side dense vector.
   TVT2  rhs_;     //!< The right-hand side dense vector.
   DRE   dres_;    //!< The dense result matrix.
   SRE   sres_;    //!< The sparse result matrix.
   ODRE  odres_;   //!< The dense result matrix with opposite storage order.
   OSRE  osres_;   //!< The sparse result matrix with opposite storage order.
   TDRE  tdres_;   //!< The transpose dense result matrix.
   TSRE  tsres_;   //!< The transpose sparse result matrix.
   TODRE todres_;  //!< The transpose dense result matrix with opposite storage order.
   TOSRE tosres_;  //!< The transpose sparse result matrix with opposite storage order.
   RT1   reflhs_;  //!< The reference left-hand side vector.
   RT2   refrhs_;  //!< The reference right-hand side vector.
   RRE   refres_;  //!< The reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT1   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT2   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TVT1  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TVT2  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( RT1   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( RT2   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RRE   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( ODRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OSRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TODRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TOSRE );

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE      ( VT1   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE      ( VT2   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE         ( TVT1  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE         ( TVT2  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE      ( RT1   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE         ( RT2   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ODRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ODRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TODRE );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TODRE );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET1, blaze::ElementType_t<TVT1>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET2, blaze::ElementType_t<TVT2>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<DRE>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<ODRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<TDRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<TODRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<SRE>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<SRE>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<OSRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<TSRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<TOSRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<DRE>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT1, blaze::TransposeType_t<TVT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT2, blaze::TransposeType_t<TVT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DRE, blaze::OppositeType_t<ODRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DRE, blaze::TransposeType_t<TDRE> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SRE, blaze::OppositeType_t<OSRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SRE, blaze::TransposeType_t<TSRE> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( OuterExprType, blaze::ResultType_t<OuterExprType>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( OuterExprType, blaze::OppositeType_t<OuterExprType>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( OuterExprType, blaze::TransposeType_t<OuterExprType> );

   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( OuterExprType, blaze::BaseType_t<OuterExprType> );
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
/*!\brief Constructor for the dense vector/dense vector outer product operation test.
//
// \param creator1 The creator for the left-hand side dense vector of the vector outer product.
// \param creator2 The creator for the right-hand side dense vector of the vector outer product.
// \exception std::runtime_error Operation error detected.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
OperationTest<VT1,VT2>::OperationTest( const Creator<VT1>& creator1, const Creator<VT2>& creator2 )
   : lhs_( creator1() )           // The left-hand side dense vector
   , rhs_( trans( creator2() ) )  // The right-hand side dense vector
   , dres_()                      // The dense result matrix
   , sres_()                      // The sparse result matrix
   , odres_()                     // The dense result matrix with opposite storage order
   , osres_()                     // The sparse result matrix with opposite storage order
   , tdres_()                     // The transpose dense result matrix
   , tsres_()                     // The transpose sparse result matrix
   , todres_()                    // The transpose dense result matrix with opposite storage order
   , tosres_()                    // The transpose sparse result matrix with opposite storage order
   , reflhs_( lhs_ )              // The reference left-hand side vector
   , refrhs_( rhs_ )              // The reference right-hand side vector
   , refres_()                    // The reference result
   , test_()                      // Label of the currently performed test
   , error_()                     // Description of the current error type
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
   testDeclSymOperation( Or_t< IsSquare<DRE>, IsResizable<DRE> >() );
   testDeclHermOperation( Or_t< IsSquare<DRE>, IsResizable<DRE> >() );
   testDeclLowOperation( And_t< Or_t< IsSquare<DRE>, IsResizable<DRE> >, Not_t< IsUniform<VT2> > >() );
   testDeclUppOperation( And_t< Or_t< IsSquare<DRE>, IsResizable<DRE> >, Not_t< IsUniform<VT1> > >() );
   testDeclDiagOperation( Or_t< IsSquare<DRE>, IsResizable<DRE> >() );
   testSubmatrixOperation( Not_t< IsUniform<DRE> >() );
   testRowOperation( Not_t< IsUniform<DRE> >() );
   testRowsOperation( Not_t< IsUniform<DRE> >() );
   testColumnOperation( Not_t< IsUniform<DRE> >() );
   testColumnsOperation( Not_t< IsUniform<DRE> >() );
   testBandOperation( Not_t< IsUniform<DRE> >() );
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
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testInitialStatus()
{
   // Checking the size of the left-hand side operand
   if( lhs_.size() != reflhs_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side dense operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
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
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Detected size = " << rhs_.size() << "\n"
          << "   Expected size = " << refrhs_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the left-hand side operand
   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of left-hand side dense operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
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
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
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
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testAssignment()
{
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
          << "   Left-hand side dense vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of left-hand side dense operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense vector type:\n"
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
          << "   Transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
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
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testEvaluation()
{
   using blaze::IsRowVector;


   {
      const auto res   ( evaluate( lhs_    * rhs_    ) );
      const auto refres( evaluate( reflhs_ * refrhs_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given vectors\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side transpose dense vector type:\n"
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
      const auto res   ( evaluate( eval( lhs_ )    * eval( rhs_ )    ) );
      const auto refres( evaluate( eval( reflhs_ ) * eval( refrhs_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated vectors\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side transpose dense vector type:\n"
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
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the matrix element access.
//
// \return void
// \exception std::runtime_error Element access error detected.
//
// This function tests the element access via the subscript operator. In case any
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testElementAccess()
{
   using blaze::equal;

   if( lhs_.size() > 0UL && rhs_.size() > 0UL )
   {
      const size_t m( lhs_.size() - 1UL );
      const size_t n( rhs_.size() - 1UL );

      if( !equal( ( lhs_ * rhs_ )(m,n), ( reflhs_ * refrhs_ )(m,n) ) ||
          !equal( ( lhs_ * rhs_ ).at(m,n), ( reflhs_ * refrhs_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of outer product expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side transpose dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( lhs_ * eval( rhs_ ) )(m,n), ( reflhs_ * eval( refrhs_ ) )(m,n) ) ||
          !equal( ( lhs_ * eval( rhs_ ) ).at(m,n), ( reflhs_ * eval( refrhs_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side transpose dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * rhs_ )(m,n), ( eval( reflhs_ ) * refrhs_ )(m,n) ) ||
          !equal( ( eval( lhs_ ) * rhs_ ).at(m,n), ( eval( reflhs_ ) * refrhs_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side transpose dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * eval( rhs_ ) )(m,n), ( eval( reflhs_ ) * eval( refrhs_ ) )(m,n) ) ||
          !equal( ( eval( lhs_ ) * eval( rhs_ ) ).at(m,n), ( eval( reflhs_ ) * eval( refrhs_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated addition expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side transpose dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      ( lhs_ * rhs_ ).at( 0UL, rhs_.size() );

      std::ostringstream oss;
      oss << " Test : Checked element access of outer product expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side dense vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      ( lhs_ * rhs_ ).at( lhs_.size(), 0UL );

      std::ostringstream oss;
      oss << " Test : Checked element access of outer product expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side dense vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the plain outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Outer product
      //=====================================================================================

      // Outer product with the given vectors
      {
         test_  = "Outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = lhs_ * rhs_;
            odres_  = lhs_ * rhs_;
            sres_   = lhs_ * rhs_;
            osres_  = lhs_ * rhs_;
            refres_ = reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Outer product with evaluated vectors
      {
         test_  = "Outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = eval( lhs_ ) * eval( rhs_ );
            odres_  = eval( lhs_ ) * eval( rhs_ );
            sres_   = eval( lhs_ ) * eval( rhs_ );
            osres_  = eval( lhs_ ) * eval( rhs_ );
            refres_ = eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Outer product with addition assignment
      //=====================================================================================

      // Outer product with addition assignment with the given vectors
      {
         test_  = "Outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += lhs_ * rhs_;
            odres_  += lhs_ * rhs_;
            sres_   += lhs_ * rhs_;
            osres_  += lhs_ * rhs_;
            refres_ += reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Outer product with addition assignment with evaluated vectors
      {
         test_  = "Outer product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += eval( lhs_ ) * eval( rhs_ );
            odres_  += eval( lhs_ ) * eval( rhs_ );
            sres_   += eval( lhs_ ) * eval( rhs_ );
            osres_  += eval( lhs_ ) * eval( rhs_ );
            refres_ += eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Outer product with subtraction assignment
      //=====================================================================================

      // Outer product with subtraction assignment with the given vectors
      {
         test_  = "Outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= lhs_ * rhs_;
            odres_  -= lhs_ * rhs_;
            sres_   -= lhs_ * rhs_;
            osres_  -= lhs_ * rhs_;
            refres_ -= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Outer product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= eval( lhs_ ) * eval( rhs_ );
            odres_  -= eval( lhs_ ) * eval( rhs_ );
            sres_   -= eval( lhs_ ) * eval( rhs_ );
            osres_  -= eval( lhs_ ) * eval( rhs_ );
            refres_ -= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Outer product with Schur product assignment
      //=====================================================================================

      // Outer product with Schur product assignment with the given vectors
      {
         test_  = "Outer product with Schur product assignment with the given vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= lhs_ * rhs_;
            odres_  %= lhs_ * rhs_;
            sres_   %= lhs_ * rhs_;
            osres_  %= lhs_ * rhs_;
            refres_ %= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Outer product with Schur product assignment with evaluated vectors
      {
         test_  = "Outer product with Schur product assignment with evaluated vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= eval( lhs_ ) * eval( rhs_ );
            odres_  %= eval( lhs_ ) * eval( rhs_ );
            sres_   %= eval( lhs_ ) * eval( rhs_ );
            osres_  %= eval( lhs_ ) * eval( rhs_ );
            refres_ %= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the negated dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the negated outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      //=====================================================================================
      // Negated outer product
      //=====================================================================================

      // Negated outer product with the given vectors
      {
         test_  = "Negated outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = -( lhs_ * rhs_ );
            odres_  = -( lhs_ * rhs_ );
            sres_   = -( lhs_ * rhs_ );
            osres_  = -( lhs_ * rhs_ );
            refres_ = -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Negated outer product with evaluated vectors
      {
         test_  = "Negated outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = -( eval( lhs_ ) * eval( rhs_ ) );
            odres_  = -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   = -( eval( lhs_ ) * eval( rhs_ ) );
            osres_  = -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Negated outer product with addition assignment
      //=====================================================================================

      // Negated outer product with addition assignment with the given vectors
      {
         test_  = "Negated outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -( lhs_ * rhs_ );
            odres_  += -( lhs_ * rhs_ );
            sres_   += -( lhs_ * rhs_ );
            osres_  += -( lhs_ * rhs_ );
            refres_ += -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Negated outer product with addition assignment with evaluated vectors
      {
         test_  = "Negated outer product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -( eval( lhs_ ) * eval( rhs_ ) );
            odres_  += -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   += -( eval( lhs_ ) * eval( rhs_ ) );
            osres_  += -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ += -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Negated outer product with subtraction assignment
      //=====================================================================================

      // Negated outer product with subtraction assignment with the given vectors
      {
         test_  = "Negated outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -( lhs_ * rhs_ );
            odres_  -= -( lhs_ * rhs_ );
            sres_   -= -( lhs_ * rhs_ );
            osres_  -= -( lhs_ * rhs_ );
            refres_ -= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Negated outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Negated outer product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -( eval( lhs_ ) * eval( rhs_ ) );
            odres_  -= -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   -= -( eval( lhs_ ) * eval( rhs_ ) );
            osres_  -= -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ -= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Negated outer product with Schur product assignment
      //=====================================================================================

      // Negated outer product with Schur product assignment with the given vectors
      {
         test_  = "Negated outer product with Schur product assignment with the given vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= -( lhs_ * rhs_ );
            odres_  %= -( lhs_ * rhs_ );
            sres_   %= -( lhs_ * rhs_ );
            osres_  %= -( lhs_ * rhs_ );
            refres_ %= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Negated outer product with Schur product assignment with evaluated vectors
      {
         test_  = "Negated outer product with Schur product assignment with evaluated vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= -( eval( lhs_ ) * eval( rhs_ ) );
            odres_  %= -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   %= -( eval( lhs_ ) * eval( rhs_ ) );
            osres_  %= -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ %= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the scaled dense vector/dense vector outer product.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the scaled outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
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
      // Self-scaling (M*=s)
      //=====================================================================================

      // Self-scaling (M*=s)
      {
         test_ = "Self-scaling (M*=s)";

         try {
            dres_   = lhs_ * rhs_;
            odres_  = dres_;
            sres_   = dres_;
            osres_  = dres_;
            refres_ = dres_;

            dres_   *= scalar;
            odres_  *= scalar;
            sres_   *= scalar;
            osres_  *= scalar;
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

         checkResults();
      }


      //=====================================================================================
      // Self-scaling (M=M*s)
      //=====================================================================================

      // Self-scaling (M=M*s)
      {
         test_ = "Self-scaling (M=M*s)";

         try {
            dres_   = lhs_ * rhs_;
            odres_  = dres_;
            sres_   = dres_;
            osres_  = dres_;
            refres_ = dres_;

            dres_   = dres_   * scalar;
            odres_  = odres_  * scalar;
            sres_   = sres_   * scalar;
            osres_  = osres_  * scalar;
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

         checkResults();
      }


      //=====================================================================================
      // Self-scaling (M=s*M)
      //=====================================================================================

      // Self-scaling (M=s*M)
      {
         test_ = "Self-scaling (M=s*M)";

         try {
            dres_   = lhs_ * rhs_;
            odres_  = dres_;
            sres_   = dres_;
            osres_  = dres_;
            refres_ = dres_;

            dres_   = scalar * dres_;
            odres_  = scalar * odres_;
            sres_   = scalar * sres_;
            osres_  = scalar * osres_;
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

         checkResults();
      }


      //=====================================================================================
      // Self-scaling (M/=s)
      //=====================================================================================

      // Self-scaling (M/=s)
      {
         test_ = "Self-scaling (M/=s)";

         try {
            dres_   = lhs_ * rhs_;
            odres_  = dres_;
            sres_   = dres_;
            osres_  = dres_;
            refres_ = dres_;

            dres_   /= scalar;
            odres_  /= scalar;
            sres_   /= scalar;
            osres_  /= scalar;
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

         checkResults();
      }


      //=====================================================================================
      // Self-scaling (M=M/s)
      //=====================================================================================

      // Self-scaling (M=M/s)
      {
         test_ = "Self-scaling (M=M/s)";

         try {
            dres_   = lhs_ * rhs_;
            odres_  = dres_;
            sres_   = dres_;
            osres_  = dres_;
            refres_ = dres_;

            dres_   = dres_   / scalar;
            odres_  = odres_  / scalar;
            sres_   = sres_   / scalar;
            osres_  = osres_  / scalar;
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

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product (s*OP)
      //=====================================================================================

      // Scaled outer product with the given vectors
      {
         test_  = "Scaled outer product with the given vectors (s*OP)";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = scalar * ( lhs_ * rhs_ );
            odres_  = scalar * ( lhs_ * rhs_ );
            sres_   = scalar * ( lhs_ * rhs_ );
            osres_  = scalar * ( lhs_ * rhs_ );
            refres_ = scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with evaluated vectors
      {
         test_  = "Scaled outer product with evaluated vectors (s*OP)";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            odres_  = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            osres_  = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product (OP*s)
      //=====================================================================================

      // Scaled outer product with the given vectors
      {
         test_  = "Scaled outer product with the given vectors (OP*s)";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = ( lhs_ * rhs_ ) * scalar;
            odres_  = ( lhs_ * rhs_ ) * scalar;
            sres_   = ( lhs_ * rhs_ ) * scalar;
            osres_  = ( lhs_ * rhs_ ) * scalar;
            refres_ = ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with evaluated vectors
      {
         test_  = "Scaled outer product with evaluated vectors (OP*s)";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            odres_  = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            osres_  = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product (OP/s)
      //=====================================================================================

      // Scaled outer product with the given vectors
      {
         test_  = "Scaled outer product with the given vectors (OP/s)";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = ( lhs_ * rhs_ ) / scalar;
            odres_  = ( lhs_ * rhs_ ) / scalar;
            sres_   = ( lhs_ * rhs_ ) / scalar;
            osres_  = ( lhs_ * rhs_ ) / scalar;
            refres_ = ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with evaluated vectors
      {
         test_  = "Scaled outer product with evaluated vectors (OP/s)";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            odres_  = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            osres_  = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with addition assignment (s*OP)
      //=====================================================================================

      // Scaled outer product with addition assignment with the given vectors
      {
         test_  = "Scaled outer product with addition assignment with the given vectors (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * ( lhs_ * rhs_ );
            odres_  += scalar * ( lhs_ * rhs_ );
            sres_   += scalar * ( lhs_ * rhs_ );
            osres_  += scalar * ( lhs_ * rhs_ );
            refres_ += scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with addition assignment with evaluated vectors
      {
         test_  = "Scaled outer product with addition assignment with evaluated vectors (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            odres_  += scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   += scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            osres_  += scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with addition assignment (OP*s)
      //=====================================================================================

      // Scaled outer product with addition assignment with the given vectors
      {
         test_  = "Scaled outer product with addition assignment with the given vectors (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( lhs_ * rhs_ ) * scalar;
            odres_  += ( lhs_ * rhs_ ) * scalar;
            sres_   += ( lhs_ * rhs_ ) * scalar;
            osres_  += ( lhs_ * rhs_ ) * scalar;
            refres_ += ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with addition assignment with evaluated vectors
      {
         test_  = "Scaled outer product with addition assignment with evaluated vectors (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            odres_  += ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   += ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            osres_  += ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with addition assignment (OP/s)
      //=====================================================================================

      // Scaled outer product with addition assignment with the given vectors
      {
         test_  = "Scaled outer product with addition assignment with the given vectors (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( lhs_ * rhs_ ) / scalar;
            odres_  += ( lhs_ * rhs_ ) / scalar;
            sres_   += ( lhs_ * rhs_ ) / scalar;
            osres_  += ( lhs_ * rhs_ ) / scalar;
            refres_ += ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with addition assignment with evaluated vectors
      {
         test_  = "Scaled outer product with addition assignment with evaluated vectors (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            odres_  += ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   += ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            osres_  += ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled outer product with subtraction assignment with the given vectors
      {
         test_  = "Scaled outer product with subtraction assignment with the given vectors (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * ( lhs_ * rhs_ );
            odres_  -= scalar * ( lhs_ * rhs_ );
            sres_   -= scalar * ( lhs_ * rhs_ );
            osres_  -= scalar * ( lhs_ * rhs_ );
            refres_ -= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Scaled outer product with subtraction assignment with evaluated vectors (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            odres_  -= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   -= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            osres_  -= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled outer product with subtraction assignment with the given vectors
      {
         test_  = "Scaled outer product with subtraction assignment with the given vectors (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( lhs_ * rhs_ ) * scalar;
            odres_  -= ( lhs_ * rhs_ ) * scalar;
            sres_   -= ( lhs_ * rhs_ ) * scalar;
            osres_  -= ( lhs_ * rhs_ ) * scalar;
            refres_ -= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Scaled outer product with subtraction assignment with evaluated vectors (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            odres_  -= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            osres_  -= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled outer product with subtraction assignment with the given vectors
      {
         test_  = "Scaled outer product with subtraction assignment with the given vectors (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( lhs_ * rhs_ ) / scalar;
            odres_  -= ( lhs_ * rhs_ ) / scalar;
            sres_   -= ( lhs_ * rhs_ ) / scalar;
            osres_  -= ( lhs_ * rhs_ ) / scalar;
            refres_ -= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Scaled outer product with subtraction assignment with evaluated vectors (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            odres_  -= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            osres_  -= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with Schur product assignment (s*OP)
      //=====================================================================================

      // Scaled outer product with Schur product assignment with the given vectors
      {
         test_  = "Scaled outer product with Schur product assignment with the given vectors (s*OP)";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= scalar * ( lhs_ * rhs_ );
            odres_  %= scalar * ( lhs_ * rhs_ );
            sres_   %= scalar * ( lhs_ * rhs_ );
            osres_  %= scalar * ( lhs_ * rhs_ );
            refres_ %= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with Schur product assignment with evaluated vectors
      {
         test_  = "Scaled outer product with Schur product assignment with evaluated vectors (s*OP)";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            odres_  %= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   %= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            osres_  %= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ %= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with Schur product assignment (OP*s)
      //=====================================================================================

      // Scaled outer product with Schur product assignment with the given vectors
      {
         test_  = "Scaled outer product with Schur product assignment with the given vectors (OP*s)";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= ( lhs_ * rhs_ ) * scalar;
            odres_  %= ( lhs_ * rhs_ ) * scalar;
            sres_   %= ( lhs_ * rhs_ ) * scalar;
            osres_  %= ( lhs_ * rhs_ ) * scalar;
            refres_ %= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with Schur product assignment with evaluated vectors
      {
         test_  = "Scaled outer product with Schur product assignment with evaluated vectors (OP*s)";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            odres_  %= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   %= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            osres_  %= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ %= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Scaled outer product with Schur product assignment (OP/s)
      //=====================================================================================

      // Scaled outer product with Schur product assignment with the given vectors
      {
         test_  = "Scaled outer product with Schur product assignment with the given vectors (OP/s)";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= ( lhs_ * rhs_ ) / scalar;
            odres_  %= ( lhs_ * rhs_ ) / scalar;
            sres_   %= ( lhs_ * rhs_ ) / scalar;
            osres_  %= ( lhs_ * rhs_ ) / scalar;
            refres_ %= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Scaled outer product with Schur product assignment with evaluated vectors
      {
         test_  = "Scaled outer product with Schur product assignment with evaluated vectors (OP/s)";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            odres_  %= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   %= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            osres_  %= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ %= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the transpose dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the transpose outer product with plain assignment. In case any
// error resulting from the outer product or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANS_OPERATION > 1 )
   {
      //=====================================================================================
      // Transpose outer product
      //=====================================================================================

      // Transpose outer product with the given vectors
      {
         test_  = "Transpose outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initTransposeResults();
            tdres_  = trans( lhs_ * rhs_ );
            todres_ = trans( lhs_ * rhs_ );
            tsres_  = trans( lhs_ * rhs_ );
            tosres_ = trans( lhs_ * rhs_ );
            refres_ = trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkTransposeResults();
      }

      // Transpose outer product with evaluated vectors
      {
         test_  = "Transpose outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initTransposeResults();
            tdres_  = trans( eval( lhs_ ) * eval( rhs_ ) );
            todres_ = trans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_  = trans( eval( lhs_ ) * eval( rhs_ ) );
            tosres_ = trans( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkTransposeResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the conjugate transpose dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the conjugate transpose outer product with plain assignment. In case
// any error resulting from the outer product or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testCTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION > 1 )
   {
      //=====================================================================================
      // Conjugate transpose outer product
      //=====================================================================================

      // Conjugate transpose outer product with the given vectors
      {
         test_  = "Conjugate transpose outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initTransposeResults();
            tdres_  = ctrans( lhs_ * rhs_ );
            todres_ = ctrans( lhs_ * rhs_ );
            tsres_  = ctrans( lhs_ * rhs_ );
            tosres_ = ctrans( lhs_ * rhs_ );
            refres_ = ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkTransposeResults();
      }

      // Conjugate transpose outer product with evaluated vectors
      {
         test_  = "Conjugate transpose outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initTransposeResults();
            tdres_  = ctrans( eval( lhs_ ) * eval( rhs_ ) );
            todres_ = ctrans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_  = ctrans( eval( lhs_ ) * eval( rhs_ ) );
            tosres_ = ctrans( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkTransposeResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the abs dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the abs outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
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
/*!\brief Testing the conjugate dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the conjugate outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
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
/*!\brief Testing the \a real dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the \a real outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
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
/*!\brief Testing the \a imag dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the \a imag outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
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
/*!\brief Testing the evaluated dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the evaluated outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
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
/*!\brief Testing the serialized dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the serialized outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
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
/*!\brief Testing the non-aliased dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the non-aliased outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
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
/*!\brief Testing the non-SIMD dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the non-SIMD outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
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
/*!\brief Testing the symmetric dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the symmetric outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testDeclSymOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION > 1 )
   {
      if( lhs_.size() != rhs_.size() )
         return;


      //=====================================================================================
      // Test-specific setup of the left-hand side operand
      //=====================================================================================

      VT1 lhs( lhs_ );
      reset( lhs );

      RT1 reflhs( lhs );


      //=====================================================================================
      // Test-specific setup of the right-hand side operand
      //=====================================================================================

      TVT2 rhs( rhs_ );
      reset( rhs );

      RT2 refrhs( rhs );


      //=====================================================================================
      // Declsym outer product
      //=====================================================================================

      // Declsym outer product with the given vectors
      {
         test_  = "Declsym outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = declsym( lhs * rhs );
            odres_  = declsym( lhs * rhs );
            sres_   = declsym( lhs * rhs );
            osres_  = declsym( lhs * rhs );
            refres_ = declsym( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Declsym outer product with evaluated vectors
      {
         test_  = "Declsym outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = declsym( eval( lhs ) * eval( rhs ) );
            odres_  = declsym( eval( lhs ) * eval( rhs ) );
            sres_   = declsym( eval( lhs ) * eval( rhs ) );
            osres_  = declsym( eval( lhs ) * eval( rhs ) );
            refres_ = declsym( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Declsym outer product with addition assignment
      //=====================================================================================

      // Declsym outer product with addition assignment with the given vectors
      {
         test_  = "Declsym outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declsym( lhs * rhs );
            odres_  += declsym( lhs * rhs );
            sres_   += declsym( lhs * rhs );
            osres_  += declsym( lhs * rhs );
            refres_ += declsym( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Declsym outer product with addition assignment with evaluated vectors
      {
         test_  = "Declsym outer product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declsym( eval( lhs ) * eval( rhs ) );
            odres_  += declsym( eval( lhs ) * eval( rhs ) );
            sres_   += declsym( eval( lhs ) * eval( rhs ) );
            osres_  += declsym( eval( lhs ) * eval( rhs ) );
            refres_ += declsym( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Declsym outer product with subtraction assignment
      //=====================================================================================

      // Declsym outer product with subtraction assignment with the given vectors
      {
         test_  = "Declsym outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declsym( lhs * rhs );
            odres_  -= declsym( lhs * rhs );
            sres_   -= declsym( lhs * rhs );
            osres_  -= declsym( lhs * rhs );
            refres_ -= declsym( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Declsym outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Declsym outer product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declsym( eval( lhs ) * eval( rhs ) );
            odres_  -= declsym( eval( lhs ) * eval( rhs ) );
            sres_   -= declsym( eval( lhs ) * eval( rhs ) );
            osres_  -= declsym( eval( lhs ) * eval( rhs ) );
            refres_ -= declsym( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Declsym outer product with Schur product assignment
      //=====================================================================================

      // Declsym outer product with Schur product assignment with the given vectors
      {
         test_  = "Declsym outer product with Schur product assignment with the given vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declsym( lhs * rhs );
            odres_  %= declsym( lhs * rhs );
            sres_   %= declsym( lhs * rhs );
            osres_  %= declsym( lhs * rhs );
            refres_ %= declsym( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Declsym outer product with Schur product assignment with evaluated vectors
      {
         test_  = "Declsym outer product with Schur product assignment with evaluated vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declsym( eval( lhs ) * eval( rhs ) );
            odres_  %= declsym( eval( lhs ) * eval( rhs ) );
            sres_   %= declsym( eval( lhs ) * eval( rhs ) );
            osres_  %= declsym( eval( lhs ) * eval( rhs ) );
            refres_ %= declsym( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the symmetric dense vector/dense vector outer product.
//
// \return void
//
// This function is called in case the symmetric vector/vector outer product operation is not
// available for the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testDeclSymOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the Hermitian dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the Hermitian outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testDeclHermOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION > 1 )
   {
      if( lhs_.size() != rhs_.size() )
         return;


      //=====================================================================================
      // Test-specific setup of the left-hand side operand
      //=====================================================================================

      VT1 lhs( lhs_ );
      reset( lhs );

      RT1 reflhs( lhs );


      //=====================================================================================
      // Test-specific setup of the right-hand side operand
      //=====================================================================================

      TVT2 rhs( rhs_ );
      reset( rhs );

      RT2 refrhs( rhs );


      //=====================================================================================
      // Declherm outer product
      //=====================================================================================

      // Declherm outer product with the given vectors
      {
         test_  = "Declherm outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = declherm( lhs * rhs );
            odres_  = declherm( lhs * rhs );
            sres_   = declherm( lhs * rhs );
            osres_  = declherm( lhs * rhs );
            refres_ = declherm( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Declherm outer product with evaluated vectors
      {
         test_  = "Declherm outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = declherm( eval( lhs ) * eval( rhs ) );
            odres_  = declherm( eval( lhs ) * eval( rhs ) );
            sres_   = declherm( eval( lhs ) * eval( rhs ) );
            osres_  = declherm( eval( lhs ) * eval( rhs ) );
            refres_ = declherm( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Declherm outer product with addition assignment
      //=====================================================================================

      // Declherm outer product with addition assignment with the given vectors
      {
         test_  = "Declherm outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declherm( lhs * rhs );
            odres_  += declherm( lhs * rhs );
            sres_   += declherm( lhs * rhs );
            osres_  += declherm( lhs * rhs );
            refres_ += declherm( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Declherm outer product with addition assignment with evaluated vectors
      {
         test_  = "Declherm outer product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declherm( eval( lhs ) * eval( rhs ) );
            odres_  += declherm( eval( lhs ) * eval( rhs ) );
            sres_   += declherm( eval( lhs ) * eval( rhs ) );
            osres_  += declherm( eval( lhs ) * eval( rhs ) );
            refres_ += declherm( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Declherm outer product with subtraction assignment
      //=====================================================================================

      // Declherm outer product with subtraction assignment with the given vectors
      {
         test_  = "Declherm outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declherm( lhs * rhs );
            odres_  -= declherm( lhs * rhs );
            sres_   -= declherm( lhs * rhs );
            osres_  -= declherm( lhs * rhs );
            refres_ -= declherm( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Declherm outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Declherm outer product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declherm( eval( lhs ) * eval( rhs ) );
            odres_  -= declherm( eval( lhs ) * eval( rhs ) );
            sres_   -= declherm( eval( lhs ) * eval( rhs ) );
            osres_  -= declherm( eval( lhs ) * eval( rhs ) );
            refres_ -= declherm( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Declherm outer product with Schur product assignment
      //=====================================================================================

      // Declherm outer product with Schur product assignment with the given vectors
      {
         test_  = "Declherm outer product with Schur product assignment with the given vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declherm( lhs * rhs );
            odres_  %= declherm( lhs * rhs );
            sres_   %= declherm( lhs * rhs );
            osres_  %= declherm( lhs * rhs );
            refres_ %= declherm( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Declherm outer product with Schur product assignment with evaluated vectors
      {
         test_  = "Declherm outer product with Schur product assignment with evaluated vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declherm( eval( lhs ) * eval( rhs ) );
            odres_  %= declherm( eval( lhs ) * eval( rhs ) );
            sres_   %= declherm( eval( lhs ) * eval( rhs ) );
            osres_  %= declherm( eval( lhs ) * eval( rhs ) );
            refres_ %= declherm( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the Hermitian dense vector/dense vector outer product.
//
// \return void
//
// This function is called in case the Hermitian vector/vector outer product operation is not
// available for the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testDeclHermOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the lower dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the lower outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testDeclLowOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLLOW_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLLOW_OPERATION > 1 )
   {
      if( lhs_.size() != rhs_.size() )
         return;


      //=====================================================================================
      // Test-specific setup of the right-hand side operand
      //=====================================================================================

      TVT2 rhs( rhs_ );
      reset( subvector( rhs, 1UL, rhs.size()-1UL ) );

      RT2 refrhs( rhs );


      //=====================================================================================
      // Decllow outer product
      //=====================================================================================

      // Decllow outer product with the given vectors
      {
         test_  = "Decllow outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = decllow( lhs_ * rhs );
            odres_  = decllow( lhs_ * rhs );
            sres_   = decllow( lhs_ * rhs );
            osres_  = decllow( lhs_ * rhs );
            refres_ = decllow( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Decllow outer product with evaluated vectors
      {
         test_  = "Decllow outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = decllow( eval( lhs_ ) * eval( rhs ) );
            odres_  = decllow( eval( lhs_ ) * eval( rhs ) );
            sres_   = decllow( eval( lhs_ ) * eval( rhs ) );
            osres_  = decllow( eval( lhs_ ) * eval( rhs ) );
            refres_ = decllow( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Decllow outer product with addition assignment
      //=====================================================================================

      // Decllow outer product with addition assignment with the given vectors
      {
         test_  = "Decllow outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += decllow( lhs_ * rhs );
            odres_  += decllow( lhs_ * rhs );
            sres_   += decllow( lhs_ * rhs );
            osres_  += decllow( lhs_ * rhs );
            refres_ += decllow( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Decllow outer product with addition assignment with evaluated vectors
      {
         test_  = "Decllow outer product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += decllow( eval( lhs_ ) * eval( rhs ) );
            odres_  += decllow( eval( lhs_ ) * eval( rhs ) );
            sres_   += decllow( eval( lhs_ ) * eval( rhs ) );
            osres_  += decllow( eval( lhs_ ) * eval( rhs ) );
            refres_ += decllow( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Decllow outer product with subtraction assignment
      //=====================================================================================

      // Decllow outer product with subtraction assignment with the given vectors
      {
         test_  = "Decllow outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= decllow( lhs_ * rhs );
            odres_  -= decllow( lhs_ * rhs );
            sres_   -= decllow( lhs_ * rhs );
            osres_  -= decllow( lhs_ * rhs );
            refres_ -= decllow( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Decllow outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Decllow outer product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= decllow( eval( lhs_ ) * eval( rhs ) );
            odres_  -= decllow( eval( lhs_ ) * eval( rhs ) );
            sres_   -= decllow( eval( lhs_ ) * eval( rhs ) );
            osres_  -= decllow( eval( lhs_ ) * eval( rhs ) );
            refres_ -= decllow( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Decllow outer product with Schur product assignment
      //=====================================================================================

      // Decllow outer product with Schur product assignment with the given vectors
      {
         test_  = "Decllow outer product with Schur product assignment with the given vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= decllow( lhs_ * rhs );
            odres_  %= decllow( lhs_ * rhs );
            sres_   %= decllow( lhs_ * rhs );
            osres_  %= decllow( lhs_ * rhs );
            refres_ %= decllow( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Decllow outer product with Schur product assignment with evaluated vectors
      {
         test_  = "Decllow outer product with Schur product assignment with evaluated vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= decllow( eval( lhs_ ) * eval( rhs ) );
            odres_  %= decllow( eval( lhs_ ) * eval( rhs ) );
            sres_   %= decllow( eval( lhs_ ) * eval( rhs ) );
            osres_  %= decllow( eval( lhs_ ) * eval( rhs ) );
            refres_ %= decllow( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the lower dense vector/dense vector outer product.
//
// \return void
//
// This function is called in case the lower vector/vector outer product operation is not
// available for the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testDeclLowOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the upper dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the upper outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testDeclUppOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLUPP_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLUPP_OPERATION > 1 )
   {
      if( lhs_.size() != rhs_.size() )
         return;


      //=====================================================================================
      // Test-specific setup of the left-hand side operand
      //=====================================================================================

      VT1 lhs( lhs_ );
      reset( subvector( lhs, 1UL, lhs.size()-1UL ) );

      RT1 reflhs( lhs );


      //=====================================================================================
      // Declupp outer product
      //=====================================================================================

      // Declupp outer product with the given vectors
      {
         test_  = "Declupp outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = declupp( lhs * rhs_ );
            odres_  = declupp( lhs * rhs_ );
            sres_   = declupp( lhs * rhs_ );
            osres_  = declupp( lhs * rhs_ );
            refres_ = declupp( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Declupp outer product with evaluated vectors
      {
         test_  = "Declupp outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = declupp( eval( lhs ) * eval( rhs_ ) );
            odres_  = declupp( eval( lhs ) * eval( rhs_ ) );
            sres_   = declupp( eval( lhs ) * eval( rhs_ ) );
            osres_  = declupp( eval( lhs ) * eval( rhs_ ) );
            refres_ = declupp( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Declupp outer product with addition assignment
      //=====================================================================================

      // Declupp outer product with addition assignment with the given vectors
      {
         test_  = "Declupp outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declupp( lhs * rhs_ );
            odres_  += declupp( lhs * rhs_ );
            sres_   += declupp( lhs * rhs_ );
            osres_  += declupp( lhs * rhs_ );
            refres_ += declupp( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Declupp outer product with addition assignment with evaluated vectors
      {
         test_  = "Declupp outer product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declupp( eval( lhs ) * eval( rhs_ ) );
            odres_  += declupp( eval( lhs ) * eval( rhs_ ) );
            sres_   += declupp( eval( lhs ) * eval( rhs_ ) );
            osres_  += declupp( eval( lhs ) * eval( rhs_ ) );
            refres_ += declupp( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Declupp outer product with subtraction assignment
      //=====================================================================================

      // Declupp outer product with subtraction assignment with the given vectors
      {
         test_  = "Declupp outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declupp( lhs * rhs_ );
            odres_  -= declupp( lhs * rhs_ );
            sres_   -= declupp( lhs * rhs_ );
            osres_  -= declupp( lhs * rhs_ );
            refres_ -= declupp( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Declupp outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Declupp outer product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declupp( eval( lhs ) * eval( rhs_ ) );
            odres_  -= declupp( eval( lhs ) * eval( rhs_ ) );
            sres_   -= declupp( eval( lhs ) * eval( rhs_ ) );
            osres_  -= declupp( eval( lhs ) * eval( rhs_ ) );
            refres_ -= declupp( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Declupp outer product with Schur product assignment
      //=====================================================================================

      // Declupp outer product with Schur product assignment with the given vectors
      {
         test_  = "Declupp outer product with Schur product assignment with the given vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declupp( lhs * rhs_ );
            odres_  %= declupp( lhs * rhs_ );
            sres_   %= declupp( lhs * rhs_ );
            osres_  %= declupp( lhs * rhs_ );
            refres_ %= declupp( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Declupp outer product with Schur product assignment with evaluated vectors
      {
         test_  = "Declupp outer product with Schur product assignment with evaluated vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declupp( eval( lhs ) * eval( rhs_ ) );
            odres_  %= declupp( eval( lhs ) * eval( rhs_ ) );
            sres_   %= declupp( eval( lhs ) * eval( rhs_ ) );
            osres_  %= declupp( eval( lhs ) * eval( rhs_ ) );
            refres_ %= declupp( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the upper dense vector/dense vector outer product.
//
// \return void
//
// This function is called in case the upper vector/vector outer product operation is not
// available for the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testDeclUppOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the diagonal dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the diagonal outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testDeclDiagOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLDIAG_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLDIAG_OPERATION > 1 )
   {
      if( lhs_.size() != rhs_.size() )
         return;


      //=====================================================================================
      // Test-specific setup of the left-hand side operand
      //=====================================================================================

      VT1 lhs( lhs_ );
      reset( lhs );

      RT1 reflhs( lhs );


      //=====================================================================================
      // Test-specific setup of the right-hand side operand
      //=====================================================================================

      TVT2 rhs( rhs_ );
      reset( rhs );

      RT2 refrhs( rhs );


      //=====================================================================================
      // Decldiag outer product
      //=====================================================================================

      // Decldiag outer product with the given vectors
      {
         test_  = "Decldiag outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = decldiag( lhs * rhs );
            odres_  = decldiag( lhs * rhs );
            sres_   = decldiag( lhs * rhs );
            osres_  = decldiag( lhs * rhs );
            refres_ = decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Decldiag outer product with evaluated vectors
      {
         test_  = "Decldiag outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            dres_   = decldiag( eval( lhs ) * eval( rhs ) );
            odres_  = decldiag( eval( lhs ) * eval( rhs ) );
            sres_   = decldiag( eval( lhs ) * eval( rhs ) );
            osres_  = decldiag( eval( lhs ) * eval( rhs ) );
            refres_ = decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Decldiag outer product with addition assignment
      //=====================================================================================

      // Decldiag outer product with addition assignment with the given vectors
      {
         test_  = "Decldiag outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += decldiag( lhs * rhs );
            odres_  += decldiag( lhs * rhs );
            sres_   += decldiag( lhs * rhs );
            osres_  += decldiag( lhs * rhs );
            refres_ += decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Decldiag outer product with addition assignment with evaluated vectors
      {
         test_  = "Decldiag outer product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += decldiag( eval( lhs ) * eval( rhs ) );
            odres_  += decldiag( eval( lhs ) * eval( rhs ) );
            sres_   += decldiag( eval( lhs ) * eval( rhs ) );
            osres_  += decldiag( eval( lhs ) * eval( rhs ) );
            refres_ += decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Decldiag outer product with subtraction assignment
      //=====================================================================================

      // Decldiag outer product with subtraction assignment with the given vectors
      {
         test_  = "Decldiag outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= decldiag( lhs * rhs );
            odres_  -= decldiag( lhs * rhs );
            sres_   -= decldiag( lhs * rhs );
            osres_  -= decldiag( lhs * rhs );
            refres_ -= decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Decldiag outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Decldiag outer product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= decldiag( eval( lhs ) * eval( rhs ) );
            odres_  -= decldiag( eval( lhs ) * eval( rhs ) );
            sres_   -= decldiag( eval( lhs ) * eval( rhs ) );
            osres_  -= decldiag( eval( lhs ) * eval( rhs ) );
            refres_ -= decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Decldiag outer product with Schur product assignment
      //=====================================================================================

      // Decldiag outer product with Schur product assignment with the given vectors
      {
         test_  = "Decldiag outer product with Schur product assignment with the given vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= decldiag( lhs * rhs );
            odres_  %= decldiag( lhs * rhs );
            sres_   %= decldiag( lhs * rhs );
            osres_  %= decldiag( lhs * rhs );
            refres_ %= decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Decldiag outer product with Schur product assignment with evaluated vectors
      {
         test_  = "Decldiag outer product with Schur product assignment with evaluated vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= decldiag( eval( lhs ) * eval( rhs ) );
            odres_  %= decldiag( eval( lhs ) * eval( rhs ) );
            sres_   %= decldiag( eval( lhs ) * eval( rhs ) );
            osres_  %= decldiag( eval( lhs ) * eval( rhs ) );
            refres_ %= decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the diagonal dense vector/dense vector outer product.
//
// \return void
//
// This function is called in case the diagonal vector/vector outer product operation is not
// available for the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testDeclDiagOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the submatrix-wise dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the submatrix-wise outer product with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the outer product or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testSubmatrixOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION > 1 )
   {
      if( lhs_.size() == 0UL || rhs_.size() == 0UL )
         return;


      //=====================================================================================
      // Submatrix-wise outer product
      //=====================================================================================

      // Submatrix-wise outer product with the given vectors
      {
         test_  = "Submatrix-wise outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( reflhs_ * refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Submatrix-wise outer product with evaluated vectors
      {
         test_  = "Submatrix-wise outer product with evaluated vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Submatrix-wise outer product with addition assignment
      //=====================================================================================

      // Submatrix-wise outer product with addition assignment with the given vectors
      {
         test_  = "Submatrix-wise outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( reflhs_ * refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Submatrix-wise outer product with addition assignment with evaluated vectors
      {
         test_  = "Submatrix-wise outer product with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Submatrix-wise outer product with subtraction assignment
      //=====================================================================================

      // Submatrix-wise outer product with subtraction assignment with the given vectors
      {
         test_  = "Submatrix-wise outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( reflhs_ * refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Submatrix-wise outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Submatrix-wise outer product with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Submatrix-wise outer product with Schur product assignment
      //=====================================================================================

      // Submatrix-wise outer product with Schur product assignment with the given vectors
      {
         test_  = "Submatrix-wise outer product with Schur product assignment with the given vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( reflhs_ * refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Submatrix-wise outer product with Schur product assignment with evaluated vectors
      {
         test_  = "Submatrix-wise outer product with Schur product assignment with evaluated vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.size() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.size() - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds submatrix construction (invalid number of rows)
      {
         test_  = "Out-of-bounds submatrix construction (invalid number of rows)";
         error_ = "Setup of out-of-bounds submatrix succeeded";

         try {
            auto sm = submatrix( lhs_ * rhs_, 1UL, 0UL, lhs_.size(), rhs_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }
      }

      // Out-of-bounds access (invalid number of columns)
      {
         test_  = "Out-of-bounds submatrix construction (invalid number of columns)";
         error_ = "Setup of out-of-bounds submatrix succeeded";

         try {
            auto sm = submatrix( lhs_ * rhs_, 0UL, 1UL, lhs_.size(), rhs_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }
      }

      // Out-of-bounds access (invalid row index)
      {
         test_  = "Out-of-bounds submatrix construction (invalid row index)";
         error_ = "Setup of out-of-bounds submatrix succeeded";

         try {
            auto sm = submatrix( lhs_ * rhs_, lhs_.size(), 0UL, 1UL, rhs_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }
      }

      // Out-of-bounds access (invalid column index)
      {
         test_  = "Out-of-bounds submatrix construction (invalid column index)";
         error_ = "Setup of out-of-bounds submatrix succeeded";

         try {
            auto sm = submatrix( lhs_ * rhs_, 0UL, rhs_.size(), lhs_.size(), 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the submatrix-wise dense vector/dense vector outer product.
//
// \return void
//
// This function is called in case the submatrix-wise outer product operation is not available
// for the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testSubmatrixOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the row-wise dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the row-wise outer product with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testRowOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_ROW_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ROW_OPERATION > 1 )
   {
      //=====================================================================================
      // Row-wise outer product
      //=====================================================================================

      // Row-wise outer product with the given vectors
      {
         test_  = "Row-wise outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) = row( lhs_ * rhs_, i );
               row( odres_ , i ) = row( lhs_ * rhs_, i );
               row( sres_  , i ) = row( lhs_ * rhs_, i );
               row( osres_ , i ) = row( lhs_ * rhs_, i );
               row( refres_, i ) = row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Row-wise outer product with evaluated vectors
      {
         test_  = "Row-wise outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) = row( eval( lhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) = row( eval( lhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) = row( eval( lhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) = row( eval( lhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) = row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Row-wise outer product with addition assignment
      //=====================================================================================

      // Row-wise outer product with addition assignment with the given vectors
      {
         test_  = "Row-wise outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) += row( lhs_ * rhs_, i );
               row( odres_ , i ) += row( lhs_ * rhs_, i );
               row( sres_  , i ) += row( lhs_ * rhs_, i );
               row( osres_ , i ) += row( lhs_ * rhs_, i );
               row( refres_, i ) += row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Row-wise outer product with addition assignment with evaluated vectors
      {
         test_  = "Row-wise outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) += row( eval( lhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) += row( eval( lhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) += row( eval( lhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) += row( eval( lhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) += row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Row-wise outer product with subtraction assignment
      //=====================================================================================

      // Row-wise outer product with subtraction assignment with the given vectors
      {
         test_  = "Row-wise outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) -= row( lhs_ * rhs_, i );
               row( odres_ , i ) -= row( lhs_ * rhs_, i );
               row( sres_  , i ) -= row( lhs_ * rhs_, i );
               row( osres_ , i ) -= row( lhs_ * rhs_, i );
               row( refres_, i ) -= row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Row-wise outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Row-wise outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) -= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) -= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) -= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) -= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) -= row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Row-wise outer product with multiplication assignment
      //=====================================================================================

      // Row-wise outer product with multiplication assignment with the given vectors
      {
         test_  = "Row-wise outer product with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) *= row( lhs_ * rhs_, i );
               row( odres_ , i ) *= row( lhs_ * rhs_, i );
               row( sres_  , i ) *= row( lhs_ * rhs_, i );
               row( osres_ , i ) *= row( lhs_ * rhs_, i );
               row( refres_, i ) *= row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Row-wise outer product with multiplication assignment with evaluated vectors
      {
         test_  = "Row-wise outer product with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.size(); ++i ) {
               row( dres_  , i ) *= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) *= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) *= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) *= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) *= row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid row index)
      {
         test_  = "Out-of-bounds row construction (invalid row index)";
         error_ = "Setup of out-of-bounds row succeeded";

         try {
            auto r = row( lhs_ * rhs_, lhs_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the row-wise dense vector/dense vector outer product.
//
// \return void
//
// This function is called in case the row-wise outer product operation is not available for the
// given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testRowOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the rows-wise dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the rows-wise outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testRowsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_ROWS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ROWS_OPERATION > 1 )
   {
      if( lhs_.size() == 0UL )
         return;


      std::vector<size_t> indices( lhs_.size() );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );


      //=====================================================================================
      // Rows-wise multiplication
      //=====================================================================================

      // Rows-wise multiplication with the given vectors
      {
         test_  = "Rows-wise multiplication with the given vectors";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( lhs_ * rhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( lhs_ * rhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( lhs_ * rhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( lhs_ * rhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Rows-wise multiplication with evaluated vectors
      {
         test_  = "Rows-wise multiplication with evaluated vectors";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Rows-wise multiplication with addition assignment
      //=====================================================================================

      // Rows-wise multiplication with addition assignment with the given vectors
      {
         test_  = "Rows-wise multiplication with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( lhs_ * rhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( lhs_ * rhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( lhs_ * rhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( lhs_ * rhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Rows-wise multiplication with addition assignment with evaluated vectors
      {
         test_  = "Rows-wise multiplication with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Rows-wise multiplication with subtraction assignment
      //=====================================================================================

      // Rows-wise multiplication with subtraction assignment with the given vectors
      {
         test_  = "Rows-wise multiplication with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( lhs_ * rhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( lhs_ * rhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( lhs_ * rhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( lhs_ * rhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Rows-wise multiplication with subtraction assignment with evaluated vectors
      {
         test_  = "Rows-wise multiplication with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Rows-wise multiplication with Schur product assignment
      //=====================================================================================

      // Rows-wise multiplication with Schur product assignment with the given vectors
      {
         test_  = "Rows-wise multiplication with Schur product assignment with the given vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( lhs_ * rhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( lhs_ * rhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( lhs_ * rhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( lhs_ * rhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Rows-wise multiplication with Schur product assignment with evaluated vectors
      {
         test_  = "Rows-wise multiplication with Schur product assignment with evaluated vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid row index; initializer_list)
      {
         test_  = "Out-of-bounds row selection construction (invalid row index; initializer_list)";
         error_ = "Setup of out-of-bounds row selection succeeded";

         try {
            auto r = rows( lhs_ * rhs_, { lhs_.size() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }
      }

      // Out-of-bounds access (invalid row index; lambda)
      {
         test_  = "Out-of-bounds row selection construction (invalid row index; lambda)";
         error_ = "Setup of out-of-bounds row selection succeeded";

         try {
            auto r = rows( lhs_ * rhs_, [index=lhs_.size()]( size_t ) { return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the rows-wise dense vector/dense vector outer product.
//
// \return void
//
// This function is called in case the rows-wise outer product operation is not available for the
// given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testRowsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the column-wise dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the column-wise outer product with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testColumnOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION
   if( BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION > 1 )
   {
      //=====================================================================================
      // Column-wise outer product
      //=====================================================================================

      // Column-wise outer product with the given vectors
      {
         test_  = "Column-wise outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) = column( lhs_ * rhs_, i );
               column( odres_ , i ) = column( lhs_ * rhs_, i );
               column( sres_  , i ) = column( lhs_ * rhs_, i );
               column( osres_ , i ) = column( lhs_ * rhs_, i );
               column( refres_, i ) = column( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Column-wise outer product with evaluated vectors
      {
         test_  = "Column-wise outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) = column( eval( lhs_ ) * eval( rhs_ ), i );
               column( odres_ , i ) = column( eval( lhs_ ) * eval( rhs_ ), i );
               column( sres_  , i ) = column( eval( lhs_ ) * eval( rhs_ ), i );
               column( osres_ , i ) = column( eval( lhs_ ) * eval( rhs_ ), i );
               column( refres_, i ) = column( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Column-wise outer product with addition assignment
      //=====================================================================================

      // Column-wise outer product with addition assignment with the given vectors
      {
         test_  = "Column-wise outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) += column( lhs_ * rhs_, i );
               column( odres_ , i ) += column( lhs_ * rhs_, i );
               column( sres_  , i ) += column( lhs_ * rhs_, i );
               column( osres_ , i ) += column( lhs_ * rhs_, i );
               column( refres_, i ) += column( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Column-wise outer product with addition assignment with evaluated vectors
      {
         test_  = "Column-wise outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) += column( eval( lhs_ ) * eval( rhs_ ), i );
               column( odres_ , i ) += column( eval( lhs_ ) * eval( rhs_ ), i );
               column( sres_  , i ) += column( eval( lhs_ ) * eval( rhs_ ), i );
               column( osres_ , i ) += column( eval( lhs_ ) * eval( rhs_ ), i );
               column( refres_, i ) += column( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Column-wise outer product with subtraction assignment
      //=====================================================================================

      // Column-wise outer product with subtraction assignment with the given vectors
      {
         test_  = "Column-wise outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) -= column( lhs_ * rhs_, i );
               column( odres_ , i ) -= column( lhs_ * rhs_, i );
               column( sres_  , i ) -= column( lhs_ * rhs_, i );
               column( osres_ , i ) -= column( lhs_ * rhs_, i );
               column( refres_, i ) -= column( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Column-wise outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Column-wise outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) -= column( eval( lhs_ ) * eval( rhs_ ), i );
               column( odres_ , i ) -= column( eval( lhs_ ) * eval( rhs_ ), i );
               column( sres_  , i ) -= column( eval( lhs_ ) * eval( rhs_ ), i );
               column( osres_ , i ) -= column( eval( lhs_ ) * eval( rhs_ ), i );
               column( refres_, i ) -= column( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Column-wise outer product with Schur product assignment
      //=====================================================================================

      // Column-wise outer product with Schur product assignment with the given vectors
      {
         test_  = "Column-wise outer product with Schur product assignment with the given vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) *= column( lhs_ * rhs_, i );
               column( odres_ , i ) *= column( lhs_ * rhs_, i );
               column( sres_  , i ) *= column( lhs_ * rhs_, i );
               column( osres_ , i ) *= column( lhs_ * rhs_, i );
               column( refres_, i ) *= column( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Column-wise outer product with Schur product assignment with evaluated vectors
      {
         test_  = "Column-wise outer product with Schur product assignment with the given vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<rhs_.size(); ++i ) {
               column( dres_  , i ) *= column( eval( lhs_ ) * eval( rhs_ ), i );
               column( odres_ , i ) *= column( eval( lhs_ ) * eval( rhs_ ), i );
               column( sres_  , i ) *= column( eval( lhs_ ) * eval( rhs_ ), i );
               column( osres_ , i ) *= column( eval( lhs_ ) * eval( rhs_ ), i );
               column( refres_, i ) *= column( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid column index)
      {
         test_  = "Out-of-bounds column construction (invalid column index)";
         error_ = "Setup of out-of-bounds column succeeded";

         try {
            auto c = column( lhs_ * rhs_, rhs_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the column-wise dense vector/dense vector outer product.
//
// \return void
//
// This function is called in case the column-wise outer product operation is not available for
// the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testColumnOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the columns-wise dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the columns-wise outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testColumnsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION > 1 )
   {
      if( rhs_.size() == 0UL )
         return;


      std::vector<size_t> indices( rhs_.size() );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );


      //=====================================================================================
      // Columns-wise multiplication
      //=====================================================================================

      // Columns-wise multiplication with the given vectors
      {
         test_  = "Columns-wise multiplication with the given vectors";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( lhs_ * rhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( lhs_ * rhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( lhs_ * rhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( lhs_ * rhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Columns-wise multiplication with evaluated vectors
      {
         test_  = "Columns-wise multiplication with evaluated vectors";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Columns-wise multiplication with addition assignment
      //=====================================================================================

      // Columns-wise multiplication with addition assignment with the given vectors
      {
         test_  = "Columns-wise multiplication with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( lhs_ * rhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( lhs_ * rhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( lhs_ * rhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( lhs_ * rhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Columns-wise multiplication with addition assignment with evaluated vectors
      {
         test_  = "Columns-wise multiplication with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Columns-wise multiplication with subtraction assignment
      //=====================================================================================

      // Columns-wise multiplication with subtraction assignment with the given vectors
      {
         test_  = "Columns-wise multiplication with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( lhs_ * rhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( lhs_ * rhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( lhs_ * rhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( lhs_ * rhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Columns-wise multiplication with subtraction assignment with evaluated vectors
      {
         test_  = "Columns-wise multiplication with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Columns-wise multiplication with Schur product assignment
      //=====================================================================================

      // Columns-wise multiplication with Schur product assignment with the given vectors
      {
         test_  = "Columns-wise multiplication with Schur product assignment with the given vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( lhs_ * rhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( lhs_ * rhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( lhs_ * rhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( lhs_ * rhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Columns-wise multiplication with Schur product assignment with evaluated vectors
      {
         test_  = "Columns-wise multiplication with Schur product assignment with evaluated vectors";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( eval( lhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid column index; initializer_list)
      {
         test_  = "Out-of-bounds column selection construction (invalid column index; initializer_list)";
         error_ = "Setup of out-of-bounds column selection succeeded";

         try {
            auto c = columns( lhs_ * rhs_, { rhs_.size() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }
      }

      // Out-of-bounds access (invalid column index; lambda)
      {
         test_  = "Out-of-bounds column selection construction (invalid column index; lambda)";
         error_ = "Setup of out-of-bounds column selection succeeded";

         try {
            auto c = columns( lhs_ * rhs_, [index=rhs_.size()]( size_t ) { return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the columns-wise dense vector/dense vector outer product.
//
// \return void
//
// This function is called in case the columns-wise outer product operation is not available for
// the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testColumnsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the band-wise dense vector/dense vector outer product.
//
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the band-wise outer product with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// outer product or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testBandOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_BAND_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BAND_OPERATION > 1 )
   {
      const ptrdiff_t ibegin( 1UL - lhs_.size() );
      const ptrdiff_t iend  ( rhs_.size() );


      //=====================================================================================
      // Band-wise outer product
      //=====================================================================================

      // Band-wise outer product with the given vectors
      {
         test_  = "Band-wise outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( lhs_ * rhs_, i );
               band( odres_ , i ) = band( lhs_ * rhs_, i );
               band( sres_  , i ) = band( lhs_ * rhs_, i );
               band( osres_ , i ) = band( lhs_ * rhs_, i );
               band( refres_, i ) = band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Band-wise outer product with evaluated vectors
      {
         test_  = "Band-wise outer product with the given vectors";
         error_ = "Failed outer product operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( eval( lhs_ ) * eval( rhs_ ), i );
               band( odres_ , i ) = band( eval( lhs_ ) * eval( rhs_ ), i );
               band( sres_  , i ) = band( eval( lhs_ ) * eval( rhs_ ), i );
               band( osres_ , i ) = band( eval( lhs_ ) * eval( rhs_ ), i );
               band( refres_, i ) = band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Band-wise outer product with addition assignment
      //=====================================================================================

      // Band-wise outer product with addition assignment with the given vectors
      {
         test_  = "Band-wise outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( lhs_ * rhs_, i );
               band( odres_ , i ) += band( lhs_ * rhs_, i );
               band( sres_  , i ) += band( lhs_ * rhs_, i );
               band( osres_ , i ) += band( lhs_ * rhs_, i );
               band( refres_, i ) += band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Band-wise outer product with addition assignment with evaluated vectors
      {
         test_  = "Band-wise outer product with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( eval( lhs_ ) * eval( rhs_ ), i );
               band( odres_ , i ) += band( eval( lhs_ ) * eval( rhs_ ), i );
               band( sres_  , i ) += band( eval( lhs_ ) * eval( rhs_ ), i );
               band( osres_ , i ) += band( eval( lhs_ ) * eval( rhs_ ), i );
               band( refres_, i ) += band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Band-wise outer product with subtraction assignment
      //=====================================================================================

      // Band-wise outer product with subtraction assignment with the given vectors
      {
         test_  = "Band-wise outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( lhs_ * rhs_, i );
               band( odres_ , i ) -= band( lhs_ * rhs_, i );
               band( sres_  , i ) -= band( lhs_ * rhs_, i );
               band( osres_ , i ) -= band( lhs_ * rhs_, i );
               band( refres_, i ) -= band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Band-wise outer product with subtraction assignment with evaluated vectors
      {
         test_  = "Band-wise outer product with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( eval( lhs_ ) * eval( rhs_ ), i );
               band( odres_ , i ) -= band( eval( lhs_ ) * eval( rhs_ ), i );
               band( sres_  , i ) -= band( eval( lhs_ ) * eval( rhs_ ), i );
               band( osres_ , i ) -= band( eval( lhs_ ) * eval( rhs_ ), i );
               band( refres_, i ) -= band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Band-wise outer product with multiplication assignment
      //=====================================================================================

      // Band-wise outer product with multiplication assignment with the given vectors
      {
         test_  = "Band-wise outer product with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( lhs_ * rhs_, i );
               band( odres_ , i ) *= band( lhs_ * rhs_, i );
               band( sres_  , i ) *= band( lhs_ * rhs_, i );
               band( osres_ , i ) *= band( lhs_ * rhs_, i );
               band( refres_, i ) *= band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }

      // Band-wise outer product with multiplication assignment with evaluated vectors
      {
         test_  = "Band-wise outer product with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( eval( lhs_ ) * eval( rhs_ ), i );
               band( odres_ , i ) *= band( eval( lhs_ ) * eval( rhs_ ), i );
               band( sres_  , i ) *= band( eval( lhs_ ) * eval( rhs_ ), i );
               band( osres_ , i ) *= band( eval( lhs_ ) * eval( rhs_ ), i );
               band( refres_, i ) *= band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException( ex );
         }

         checkResults();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid lower band index)
      {
         test_  = "Out-of-bounds band construction (invalid lower band index)";
         error_ = "Setup of out-of-bounds band succeeded";

         try {
            auto b = band( lhs_ * rhs_, -( lhs_.size() == 0UL ? 1UL : lhs_.size() ) );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }
      }

      // Out-of-bounds access (invalid upper band index)
      {
         test_  = "Out-of-bounds band construction (invalid upper band index)";
         error_ = "Setup of out-of-bounds band succeeded";

         try {
            auto b = band( lhs_ * rhs_, ( rhs_.size() == 0UL ? 1UL : rhs_.size() ) );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense vector type:\n"
                << "     " << typeid( VT1 ).name() << "\n"
                << "   Right-hand side dense vector type:\n"
                << "     " << typeid( TVT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the band-wise dense vector/dense vector outer product.
//
// \return void
//
// This function is called in case the band-wise outer product operation is not available for
// the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testBandOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the customized dense vector/dense vector outer product.
//
// \param op The custom operation to be tested.
// \param name The human-readable name of the operation.
// \return void
// \exception std::runtime_error Outer product error detected.
//
// This function tests the outer product with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment in combination with a custom operation.
// In case any error resulting from the outer product or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
template< typename OP >   // Type of the custom operation
void OperationTest<VT1,VT2>::testCustomOperation( OP op, const std::string& name )
{
   //=====================================================================================
   // Customized outer product
   //=====================================================================================

   // Customized outer product with the given vectors
   {
      test_  = "Customized outer product with the given vectors (" + name + ")";
      error_ = "Failed outer product operation";

      try {
         initResults();
         dres_   = op( lhs_ * rhs_ );
         odres_  = op( lhs_ * rhs_ );
         sres_   = op( lhs_ * rhs_ );
         osres_  = op( lhs_ * rhs_ );
         refres_ = op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException( ex );
      }

      checkResults();
   }

   // Customized outer product with evaluated vectors
   {
      test_  = "Customized outer product with evaluated vectors (" + name + ")";
      error_ = "Failed outer product operation";

      try {
         initResults();
         dres_   = op( eval( lhs_ ) * eval( rhs_ ) );
         odres_  = op( eval( lhs_ ) * eval( rhs_ ) );
         sres_   = op( eval( lhs_ ) * eval( rhs_ ) );
         osres_  = op( eval( lhs_ ) * eval( rhs_ ) );
         refres_ = op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException( ex );
      }

      checkResults();
   }


   //=====================================================================================
   // Customized outer product with addition assignment
   //=====================================================================================

   // Customized outer product with addition assignment with the given vectors
   {
      test_  = "Customized outer product with addition assignment with the given vectors (" + name + ")";
      error_ = "Failed addition assignment operation";

      try {
         initResults();
         dres_   += op( lhs_ * rhs_ );
         odres_  += op( lhs_ * rhs_ );
         sres_   += op( lhs_ * rhs_ );
         osres_  += op( lhs_ * rhs_ );
         refres_ += op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException( ex );
      }

      checkResults();
   }

   // Customized outer product with addition assignment with evaluated vectors
   {
      test_  = "Customized outer product with addition assignment with evaluated vectors (" + name + ")";
      error_ = "Failed addition assignment operation";

      try {
         initResults();
         dres_   += op( eval( lhs_ ) * eval( rhs_ ) );
         odres_  += op( eval( lhs_ ) * eval( rhs_ ) );
         sres_   += op( eval( lhs_ ) * eval( rhs_ ) );
         osres_  += op( eval( lhs_ ) * eval( rhs_ ) );
         refres_ += op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException( ex );
      }

      checkResults();
   }


   //=====================================================================================
   // Customized outer product with subtraction assignment
   //=====================================================================================

   // Customized outer product with subtraction assignment with the given vectors
   {
      test_  = "Customized outer product with subtraction assignment with the given vectors (" + name + ")";
      error_ = "Failed subtraction assignment operation";

      try {
         initResults();
         dres_   -= op( lhs_ * rhs_ );
         odres_  -= op( lhs_ * rhs_ );
         sres_   -= op( lhs_ * rhs_ );
         osres_  -= op( lhs_ * rhs_ );
         refres_ -= op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException( ex );
      }

      checkResults();
   }

   // Customized outer product with subtraction assignment with evaluated vectors
   {
      test_  = "Customized outer product with subtraction assignment with evaluated vectors (" + name + ")";
      error_ = "Failed subtraction assignment operation";

      try {
         initResults();
         dres_   -= op( eval( lhs_ ) * eval( rhs_ ) );
         odres_  -= op( eval( lhs_ ) * eval( rhs_ ) );
         sres_   -= op( eval( lhs_ ) * eval( rhs_ ) );
         osres_  -= op( eval( lhs_ ) * eval( rhs_ ) );
         refres_ -= op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException( ex );
      }

      checkResults();
   }


   //=====================================================================================
   // Customized outer product with Schur product assignment
   //=====================================================================================

   // Customized outer product with Schur product assignment with the given vectors
   {
      test_  = "Customized outer product with Schur product assignment with the given vectors (" + name + ")";
      error_ = "Failed Schur product assignment operation";

      try {
         initResults();
         dres_   %= op( lhs_ * rhs_ );
         odres_  %= op( lhs_ * rhs_ );
         sres_   %= op( lhs_ * rhs_ );
         osres_  %= op( lhs_ * rhs_ );
         refres_ %= op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException( ex );
      }

      checkResults();
   }

   // Customized outer product with Schur product assignment with evaluated vectors
   {
      test_  = "Customized outer product with Schur product assignment with evaluated vectors (" + name + ")";
      error_ = "Failed Schur product assignment operation";

      try {
         initResults();
         dres_   %= op( eval( lhs_ ) * eval( rhs_ ) );
         odres_  %= op( eval( lhs_ ) * eval( rhs_ ) );
         sres_   %= op( eval( lhs_ ) * eval( rhs_ ) );
         osres_  %= op( eval( lhs_ ) * eval( rhs_ ) );
         refres_ %= op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException( ex );
      }

      checkResults();
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
// This function is called after each test case to check and compare the computed results.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::checkResults()
{
   if( !isEqual( dres_, refres_ ) || !isEqual( odres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side dense vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Result:\n" << dres_ << "\n"
          << "   Result with opposite storage order:\n" << odres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( sres_, refres_ ) || !isEqual( osres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side dense vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Result:\n" << sres_ << "\n"
          << "   Result with opposite storage order:\n" << osres_ << "\n"
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
// results.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::checkTransposeResults()
{
   if( !isEqual( tdres_, refres_ ) || !isEqual( todres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect transpose dense result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side dense vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Transpose result:\n" << tdres_ << "\n"
          << "   Transpose result with opposite storage order:\n" << todres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( tsres_, refres_ ) || !isEqual( tosres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect transpose sparse result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side dense vector type:\n"
          << "     " << typeid( VT1 ).name() << "\n"
          << "   Right-hand side transpose dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n"
          << "   Transpose result:\n" << tsres_ << "\n"
          << "   Transpose result with opposite storage order:\n" << tosres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking and comparing the error message of the given exception.
//
// \param ex The exception to be checked.
// \param message The expected error message.
// \return void
// \exception std::runtime_error Wrong error message.
//
// This function is called to check the error message of the given exception. In case the error
// message does not correspond to the expected message, a \a std::runtime_error  exception is
// thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::checkExceptionMessage( const std::exception& ex, const std::string& message )
{
   if( ex.what() != message ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Wrong error message\n"
          << " Details:\n"
          << "   Error message: \"" << ex.what() << "\"\n"
          << "   Expected error message: \"" << message << "\"\n";
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
/*!\brief Initializing the non-transpose result matrices.
//
// \return void
//
// This function is called before each non-transpose test case to initialize the according result
// matrices to random values.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::initResults()
{
   const blaze::UnderlyingBuiltin_t<DRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<DRE> max( randmax );

   resize( dres_, size( lhs_ ), size( rhs_ ) );
   randomize( dres_, min, max );

   odres_  = dres_;
   sres_   = dres_;
   osres_  = dres_;
   refres_ = dres_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initializing the transpose result matrices.
//
// \return void
//
// This function is called before each transpose test case to initialize the according result
// matrices to random values.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::initTransposeResults()
{
   const blaze::UnderlyingBuiltin_t<TDRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<TDRE> max( randmax );

   resize( tdres_, size( rhs_ ), size( lhs_ ) );
   randomize( tdres_, min, max );

   todres_ = tdres_;
   tsres_  = tdres_;
   tosres_ = tdres_;
   refres_ = tdres_;
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
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::convertException( const std::exception& ex )
{
   std::ostringstream oss;
   oss << " Test : " << test_ << "\n"
       << " Error: " << error_ << "\n"
       << " Details:\n"
       << "   Random seed = " << blaze::getSeed() << "\n"
       << "   Left-hand side dense vector type:\n"
       << "     " << typeid( VT1 ).name() << "\n"
       << "   Right-hand side transpose dense vector type:\n"
       << "     " << typeid( TVT2 ).name() << "\n"
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
/*!\brief Testing the vector outer product between two specific vector types.
//
// \param creator1 The creator for the left-hand side vector.
// \param creator2 The creator for the right-hand side vector.
// \return void
*/
template< typename VT1    // Type of the left-hand side dense vector
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
//  MACROS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the definition of a dense vector/dense vector outer product test case.
*/
#define DEFINE_DVECDVECOUTER_OPERATION_TEST( VT1, VT2 ) \
   extern template class blazetest::mathtest::operations::dvecdvecouter::OperationTest<VT1,VT2>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a dense vector/dense vector outer product test case.
*/
#define RUN_DVECDVECOUTER_OPERATION_TEST( C1, C2 ) \
   blazetest::mathtest::operations::dvecdvecouter::runTest( C1, C2 )
/*! \endcond */
//*************************************************************************************************

} // namespace dvecdvecouter

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
