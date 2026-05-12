//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/dmatsvecmult/OperationTest.h
//  \brief Header file for the dense matrix/sparse vector multiplication operation test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_DMATSVECMULT_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_DMATSVECMULT_OPERATIONTEST_H_


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
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/Functors.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDivisor.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
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

namespace dmatsvecmult {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the dense matrix/sparse vector multiplication operation test.
//
// This class template represents one particular matrix/vector multiplication test between a
// matrix and a vector of particular types. The two template arguments \a MT and \a VT represent
// the types of the left-hand side matrix and right-hand side vector, respectively.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
class OperationTest
{
 private:
   //**Type definitions****************************************************************************
   using MET = blaze::ElementType_t<MT>;  //!< Element type of the matrix type
   using VET = blaze::ElementType_t<VT>;  //!< Element type of the vector type

   using OMT  = blaze::OppositeType_t<MT>;    //!< Matrix type with opposite storage order
   using TMT  = blaze::TransposeType_t<MT>;   //!< Transpose matrix type
   using TOMT = blaze::TransposeType_t<OMT>;  //!< Transpose matrix type with opposite storage order
   using TVT  = blaze::TransposeType_t<VT>;   //!< Transpose vector type

   //! Dense result type
   using DRE = blaze::MultTrait_t<MT,VT>;

   using DET  = blaze::ElementType_t<DRE>;    //!< Element type of the dense result
   using TDRE = blaze::TransposeType_t<DRE>;  //!< Transpose dense result type

   //! Sparse result type
   using SRE = blaze::CompressedVector<DET,false>;

   using SET  = blaze::ElementType_t<SRE>;    //!< Element type of the sparse result
   using TSRE = blaze::TransposeType_t<SRE>;  //!< Transpose sparse result type

   using MRT  = blaze::DynamicMatrix<MET,false>;  //!< Matrix reference type
   using VRT  = blaze::DynamicVector<VET,false>;  //!< Vector reference type
   using RRE  = blaze::MultTrait_t<MRT,VRT>;      //!< Reference result type
   using TRRE = blaze::TransposeType_t<RRE>;      //!< Transpose reference result type

   //! Type of the matrix/vector multiplication expression
   using MatVecMultExprType =
      blaze::RemoveCVRef_t< decltype( std::declval<MT>() * std::declval<VT>() ) >;

   //! Type of the transpose matrix/vector multiplication expression
   using TMatVecMultExprType =
      blaze::RemoveCVRef_t< decltype( std::declval<OMT>() * std::declval<VT>() ) >;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit OperationTest( const Creator<MT>& creator1, const Creator<VT>& creator2 );
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
   template< typename LT > void checkResults();
   template< typename LT > void checkTransposeResults();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initResults();
   void initTransposeResults();
   template< typename LT > void convertException( const std::exception& ex );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MT   lhs_;      //!< The left-hand side dense matrix.
   VT   rhs_;      //!< The right-hand side sparse vector.
   DRE  dres_;     //!< The dense result vector.
   SRE  sres_;     //!< The sparse result vector.
   MRT  reflhs_;   //!< The reference left-hand side matrix.
   VRT  refrhs_;   //!< The reference right-hand side vector.
   RRE  refres_;   //!< The reference result.
   OMT  olhs_;     //!< The left-hand side dense matrix with opposite storage order.
   TDRE tdres_;    //!< The transpose dense result vector.
   TSRE tsres_;    //!< The transpose sparse result vector.
   TRRE trefres_;  //!< The transpose reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MT   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TMT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TOMT );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TVT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MRT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VRT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( RRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( DRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( SRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TDRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TSRE );

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( MT   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TMT  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TOMT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE      ( VT   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE         ( TVT  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( MRT  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE      ( VRT  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE      ( DRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE      ( SRE  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE         ( TDRE );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE         ( TSRE );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MET, blaze::ElementType_t<OMT>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MET, blaze::ElementType_t<TMT>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MET, blaze::ElementType_t<TOMT>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VET, blaze::ElementType_t<TVT>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<DRE>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<TDRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<SRE>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<SRE>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<TSRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<DRE>    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT , blaze::OppositeType_t<OMT>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT , blaze::TransposeType_t<TMT>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT , blaze::TransposeType_t<TVT>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DRE, blaze::TransposeType_t<TDRE> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SRE, blaze::TransposeType_t<TSRE> );

   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG     ( MatVecMultExprType, blaze::ResultType_t<MatVecMultExprType>    );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_DIFFERENT_TRANSPOSE_FLAG( MatVecMultExprType, blaze::TransposeType_t<MatVecMultExprType> );

   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG     ( TMatVecMultExprType, blaze::ResultType_t<TMatVecMultExprType>    );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_DIFFERENT_TRANSPOSE_FLAG( TMatVecMultExprType, blaze::TransposeType_t<TMatVecMultExprType> );

   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( MatVecMultExprType , blaze::BaseType_t<MatVecMultExprType > );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TMatVecMultExprType, blaze::BaseType_t<TMatVecMultExprType> );
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
/*!\brief Constructor for the dense matrix/sparse vector multiplication operation test.
//
// \param creator1 The creator for the left-hand side dense matrix of the multiplication.
// \param creator2 The creator for the right-hand side sparse vector of the multiplication.
// \exception std::runtime_error Operation error detected.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
OperationTest<MT,VT>::OperationTest( const Creator<MT>& creator1, const Creator<VT>& creator2 )
   : lhs_ ( creator1() )  // The left-hand side dense matrix
   , rhs_ ( creator2() )  // The right-hand side sparse vector
   , dres_()              // The dense result vector
   , sres_()              // The sparse result vector
   , reflhs_( lhs_ )      // The reference left-hand side matrix
   , refrhs_( rhs_ )      // The reference right-hand side vector
   , refres_()            // The reference result
   , olhs_( lhs_ )        // The left-hand side dense matrix with opposite storage order.
   , tdres_()             // The transpose dense result vector.
   , tsres_()             // The transpose sparse result vector.
   , trefres_()           // The transpose reference result.
   , test_()              // Label of the currently performed test
   , error_()             // Description of the current error type
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
/*!\brief Tests on the initial status of the operands.
//
// \return void
// \exception std::runtime_error Initialization error detected.
//
// This function runs tests on the initial status of the operands. In case any initialization
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testInitialStatus()
{
   //=====================================================================================
   // Performing initial tests with the given types
   //=====================================================================================

   // Checking the number of rows of the left-hand side operand
   if( lhs_.rows() != reflhs_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side dense operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Detected number of rows = " << lhs_.rows() << "\n"
          << "   Expected number of rows = " << reflhs_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the left-hand side operand
   if( lhs_.columns() != reflhs_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side dense operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Detected number of columns = " << lhs_.columns() << "\n"
          << "   Expected number of columns = " << reflhs_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the size of the right-hand side operand
   if( rhs_.size() != refrhs_.size() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side sparse operand\n"
          << " Error: Invalid vector size\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Detected size = " << rhs_.size() << "\n"
          << "   Expected size = " << refrhs_.size() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the left-hand side operand
   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of left-hand side dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Current initialization:\n" << lhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the right-hand side operand
   if( !isEqual( rhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of right-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing initial tests with the transpose types
   //=====================================================================================

   // Checking the number of rows of the transpose left-hand side operand
   if( olhs_.rows() != reflhs_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of transpose left-hand side dense operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose dense matrix type:\n"
          << "     " << typeid( TMT ).name() << "\n"
          << "   Detected number of rows = " << olhs_.rows() << "\n"
          << "   Expected number of rows = " << reflhs_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the transpose left-hand side operand
   if( olhs_.columns() != reflhs_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of transpose left-hand side dense operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose dense matrix type:\n"
          << "     " << typeid( TMT ).name() << "\n"
          << "   Detected number of columns = " << olhs_.columns() << "\n"
          << "   Expected number of columns = " << reflhs_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the transpose left-hand side operand
   if( !isEqual( olhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of transpose left-hand side dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose dense matrix type:\n"
          << "     " << typeid( TMT ).name() << "\n"
          << "   Current initialization:\n" << olhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
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
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testAssignment()
{
   //=====================================================================================
   // Performing an assignment with the given types
   //=====================================================================================

   try {
      lhs_ = reflhs_;
      rhs_ = refrhs_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the given types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side dense matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of left-hand side dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Current initialization:\n" << lhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( rhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of right-hand side sparse operand\n"
          << " Error: Invalid vector initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing an assignment with the transpose types
   //=====================================================================================

   try {
      olhs_ = reflhs_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the transpose types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose left-hand side dense matrix type:\n"
          << "     " << typeid( TMT ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( olhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of transpose left-hand side dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose dense matrix type:\n"
          << "     " << typeid( TMT ).name() << "\n"
          << "   Current initialization:\n" << olhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
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
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testEvaluation()
{
   using blaze::IsRowMajorMatrix;


   //=====================================================================================
   // Testing the evaluation with the given types
   //=====================================================================================

   {
      const auto res   ( evaluate( lhs_    * rhs_    ) );
      const auto refres( evaluate( reflhs_ * refrhs_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrix/vector\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<MT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
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
         oss << " Test: Evaluation with evaluated matrix/vector\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<MT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
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
      const auto res   ( evaluate( olhs_   * rhs_    ) );
      const auto refres( evaluate( reflhs_ * refrhs_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the transpose matrix/vector\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<OMT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( olhs_ ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
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
      const auto res   ( evaluate( eval( olhs_ )   * eval( rhs_ )    ) );
      const auto refres( evaluate( eval( reflhs_ ) * eval( refrhs_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated transpose matrix/vector\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<OMT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( olhs_ ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
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
/*!\brief Testing the vector element access.
//
// \return void
// \exception std::runtime_error Element access error detected.
//
// This function tests the element access via the subscript operator. In case any
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testElementAccess()
{
   using blaze::equal;


   //=====================================================================================
   // Testing the element access with the given types
   //=====================================================================================

   if( lhs_.rows() > 0UL )
   {
      const size_t n( lhs_.rows() - 1UL );

      if( !equal( ( lhs_ * rhs_ )[n], ( reflhs_ * refrhs_ )[n] ) ||
          !equal( ( lhs_ * rhs_ ).at(n), ( reflhs_ * refrhs_ ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of multiplication expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( lhs_ * eval( rhs_ ) )[n], ( reflhs_ * eval( refrhs_ ) )[n] ) ||
          !equal( ( lhs_ * eval( rhs_ ) ).at(n), ( reflhs_ * eval( refrhs_ ) ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * rhs_ )[n], ( eval( reflhs_ ) * refrhs_ )[n] ) ||
          !equal( ( eval( lhs_ ) * rhs_ ).at(n), ( eval( reflhs_ ) * refrhs_ ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * eval( rhs_ ) )[n], ( eval( reflhs_ ) * eval( refrhs_ ) )[n] ) ||
          !equal( ( eval( lhs_ ) * eval( rhs_ ) ).at(n), ( eval( reflhs_ ) * eval( refrhs_ ) ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      ( lhs_ * rhs_ ).at( lhs_.rows() );

      std::ostringstream oss;
      oss << " Test : Checked element access of multiplication expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side row-major dense matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}


   //=====================================================================================
   // Testing the element access with the transpose types
   //=====================================================================================

   if( olhs_.rows() > 0UL )
   {
      const size_t n( olhs_.rows() - 1UL );

      if( !equal( ( olhs_ * rhs_ )[n], ( reflhs_ * refrhs_ )[n] ) ||
          !equal( ( olhs_ * rhs_ ).at(n), ( reflhs_ * refrhs_ ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of transpose multiplication expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( TMT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( olhs_ * eval( rhs_ ) )[n], ( reflhs_ * eval( refrhs_ ) )[n] ) ||
          !equal( ( olhs_ * eval( rhs_ ) ).at(n), ( reflhs_ * eval( refrhs_ ) ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated transpose multiplication expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( TMT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( olhs_ ) * rhs_ )[n], ( eval( reflhs_ ) * refrhs_ )[n] ) ||
          !equal( ( eval( olhs_ ) * rhs_ ).at(n), ( eval( reflhs_ ) * refrhs_ ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated transpose multiplication expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( TMT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( olhs_ ) * eval( rhs_ ) )[n], ( eval( reflhs_ ) * eval( refrhs_ ) )[n] ) ||
          !equal( ( eval( olhs_ ) * eval( rhs_ ) ).at(n), ( eval( reflhs_ ) * eval( refrhs_ ) ).at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated transpose multiplication expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( TMT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      ( olhs_ * rhs_ ).at( olhs_.rows() );

      std::ostringstream oss;
      oss << " Test : Checked element access of transpose multiplication expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side column-major dense matrix type:\n"
          << "     " << typeid( TMT ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the plain matrix/vector multiplication with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the multiplication or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Multiplication
      //=====================================================================================

      // Multiplication with the given matrix/vector
      {
         test_  = "Multiplication with the given matrix/vector";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = lhs_ * rhs_;
            sres_   = lhs_ * rhs_;
            refres_ = reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = olhs_ * rhs_;
            sres_   = olhs_ * rhs_;
            refres_ = reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Multiplication with evaluated matrix/vector
      {
         test_  = "Multiplication with evaluated matrix/vector";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = eval( lhs_ ) * eval( rhs_ );
            sres_   = eval( lhs_ ) * eval( rhs_ );
            refres_ = eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = eval( olhs_ ) * eval( rhs_ );
            sres_   = eval( olhs_ ) * eval( rhs_ );
            refres_ = eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Multiplication with addition assignment
      //=====================================================================================

      // Multiplication with addition assignment with the given matrix/vector
      {
         test_  = "Multiplication with addition assignment with the given matrix/vector";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += lhs_ * rhs_;
            sres_   += lhs_ * rhs_;
            refres_ += reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += olhs_ * rhs_;
            sres_   += olhs_ * rhs_;
            refres_ += reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Multiplication with addition assignment with evaluated matrix/vector
      {
         test_  = "Multiplication with addition assignment with evaluated matrix/vector";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += eval( lhs_ ) * eval( rhs_ );
            sres_   += eval( lhs_ ) * eval( rhs_ );
            refres_ += eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += eval( olhs_ ) * eval( rhs_ );
            sres_   += eval( olhs_ ) * eval( rhs_ );
            refres_ += eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Multiplication with subtraction assignment
      //=====================================================================================

      // Multiplication with subtraction assignment with the given matrix/vector
      {
         test_  = "Multiplication with subtraction assignment with the given matrix/vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= lhs_ * rhs_;
            sres_   -= lhs_ * rhs_;
            refres_ -= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= olhs_ * rhs_;
            sres_   -= olhs_ * rhs_;
            refres_ -= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Multiplication with subtraction assignment with evaluated matrix/vector
      {
         test_  = "Multiplication with subtraction assignment with evaluated matrix/vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= eval( lhs_ ) * eval( rhs_ );
            sres_   -= eval( lhs_ ) * eval( rhs_ );
            refres_ -= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= eval( olhs_ ) * eval( rhs_ );
            sres_   -= eval( olhs_ ) * eval( rhs_ );
            refres_ -= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Multiplication with multiplication assignment
      //=====================================================================================

      // Multiplication with multiplication assignment with the given matrix/vector
      {
         test_  = "Multiplication with multiplication assignment with the given matrix/vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= lhs_ * rhs_;
            sres_   *= lhs_ * rhs_;
            refres_ *= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   *= olhs_ * rhs_;
            sres_   *= olhs_ * rhs_;
            refres_ *= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Multiplication with multiplication assignment with evaluated matrix/vector
      {
         test_  = "Multiplication with multiplication assignment with evaluated matrix/vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= eval( lhs_ ) * eval( rhs_ );
            sres_   *= eval( lhs_ ) * eval( rhs_ );
            refres_ *= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   *= eval( olhs_ ) * eval( rhs_ );
            sres_   *= eval( olhs_ ) * eval( rhs_ );
            refres_ *= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Multiplication with division assignment
      //=====================================================================================

      if( !blaze::IsUniform_v<MT> && blaze::isDivisor( lhs_ * rhs_ ) )
      {
         // Multiplication with division assignment with the given matrix/vector
         {
            test_  = "Multiplication with division assignment with the given matrix/vector";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= lhs_ * rhs_;
               sres_   /= lhs_ * rhs_;
               refres_ /= reflhs_ * refrhs_;
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               dres_   /= olhs_ * rhs_;
               sres_   /= olhs_ * rhs_;
               refres_ /= reflhs_ * refrhs_;
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }

         // Multiplication with division assignment with evaluated matrix/vector
         {
            test_  = "Multiplication with division assignment with evaluated matrix/vector";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= eval( lhs_ ) * eval( rhs_ );
               sres_   /= eval( lhs_ ) * eval( rhs_ );
               refres_ /= eval( reflhs_ ) * eval( refrhs_ );
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               dres_   /= eval( olhs_ ) * eval( rhs_ );
               sres_   /= eval( olhs_ ) * eval( rhs_ );
               refres_ /= eval( reflhs_ ) * eval( refrhs_ );
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the negated dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the negated matrix/vector multiplication with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the multiplication or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      //=====================================================================================
      // Negated multiplication
      //=====================================================================================

      // Negated multiplication with the given matrix/vector
      {
         test_  = "Negated multiplication with the given matrix/vector";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = -( lhs_ * rhs_ );
            sres_   = -( lhs_ * rhs_ );
            refres_ = -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = -( olhs_ * rhs_ );
            sres_   = -( olhs_ * rhs_ );
            refres_ = -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Negated multiplication with evaluated matrix/vector
      {
         test_  = "Negated multiplication with evaluated matrix/vector";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   = -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = -( eval( olhs_ ) * eval( rhs_ ) );
            sres_   = -( eval( olhs_ ) * eval( rhs_ ) );
            refres_ = -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Negated multiplication with addition assignment
      //=====================================================================================

      // Negated multiplication with addition assignment with the given matrix/vector
      {
         test_  = "Negated multiplication with addition assignment with the given matrix/vector";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -( lhs_ * rhs_ );
            sres_   += -( lhs_ * rhs_ );
            refres_ += -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += -( olhs_ * rhs_ );
            sres_   += -( olhs_ * rhs_ );
            refres_ += -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Negated multiplication with addition assignment with evaluated matrix/vector
      {
         test_  = "Negated multiplication with addition assignment with evaluated matrix/vector";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   += -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ += -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += -( eval( olhs_ ) * eval( rhs_ ) );
            sres_   += -( eval( olhs_ ) * eval( rhs_ ) );
            refres_ += -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Negated multiplication with subtraction assignment
      //=====================================================================================

      // Negated multiplication with subtraction assignment with the given matrix/vector
      {
         test_  = "Negated multiplication with subtraction assignment with the given matrix/vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -( lhs_ * rhs_ );
            sres_   -= -( lhs_ * rhs_ );
            refres_ -= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= -( olhs_ * rhs_ );
            sres_   -= -( olhs_ * rhs_ );
            refres_ -= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Negated multiplication with subtraction assignment with evaluated matrix/vector
      {
         test_  = "Negated multiplication with subtraction assignment with evaluated matrix/vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   -= -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ -= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= -( eval( olhs_ ) * eval( rhs_ ) );
            sres_   -= -( eval( olhs_ ) * eval( rhs_ ) );
            refres_ -= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Negated multiplication with multiplication assignment
      //=====================================================================================

      // Negated multiplication with multiplication assignment with the given matrix/vector
      {
         test_  = "Negated multiplication with multiplication assignment with the given matrix/vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= -( lhs_ * rhs_ );
            sres_   *= -( lhs_ * rhs_ );
            refres_ *= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   *= -( olhs_ * rhs_ );
            sres_   *= -( olhs_ * rhs_ );
            refres_ *= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Negated multiplication with multiplication assignment with evaluated matrix/vector
      {
         test_  = "Negated multiplication with multiplication assignment with evaluated matrix/vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   *= -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ *= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   *= -( eval( olhs_ ) * eval( rhs_ ) );
            sres_   *= -( eval( olhs_ ) * eval( rhs_ ) );
            refres_ *= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Negated multiplication with division assignment
      //=====================================================================================

      if( !blaze::IsUniform_v<MT> && blaze::isDivisor( lhs_ * rhs_ ) )
      {
         // Negated multiplication with division assignment with the given matrix/vector
         {
            test_  = "Negated multiplication with division assignment with the given matrix/vector";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= -( lhs_ * rhs_ );
               sres_   /= -( lhs_ * rhs_ );
               refres_ /= -( reflhs_ * refrhs_ );
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               dres_   /= -( olhs_ * rhs_ );
               sres_   /= -( olhs_ * rhs_ );
               refres_ /= -( reflhs_ * refrhs_ );
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }

         // Negated multiplication with division assignment with evaluated matrix/vector
         {
            test_  = "Negated multiplication with division assignment with evaluated matrix/vector";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= -( eval( lhs_ ) * eval( rhs_ ) );
               sres_   /= -( eval( lhs_ ) * eval( rhs_ ) );
               refres_ /= -( eval( reflhs_ ) * eval( refrhs_ ) );
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               dres_   /= -( eval( olhs_ ) * eval( rhs_ ) );
               sres_   /= -( eval( olhs_ ) * eval( rhs_ ) );
               refres_ /= -( eval( reflhs_ ) * eval( refrhs_ ) );
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the scaled dense matrix/sparse vector multiplication.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the scaled matrix/vector multiplication with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the multiplication or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
template< typename T >   // Type of the scalar
void OperationTest<MT,VT>::testScaledOperation( T scalar )
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
            dres_   = lhs_ * rhs_;
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

         checkResults<MT>();
      }


      //=====================================================================================
      // Self-scaling (v=v*s)
      //=====================================================================================

      // Self-scaling (v=v*s)
      {
         test_ = "Self-scaling (v=v*s)";

         try {
            dres_   = lhs_ * rhs_;
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

         checkResults<MT>();
      }


      //=====================================================================================
      // Self-scaling (v=s*v)
      //=====================================================================================

      // Self-scaling (v=s*v)
      {
         test_ = "Self-scaling (v=s*v)";

         try {
            dres_   = lhs_ * rhs_;
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

         checkResults<MT>();
      }


      //=====================================================================================
      // Self-scaling (v/=s)
      //=====================================================================================

      // Self-scaling (v/=s)
      {
         test_ = "Self-scaling (v/=s)";

         try {
            dres_   = lhs_ * rhs_;
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

         checkResults<MT>();
      }


      //=====================================================================================
      // Self-scaling (v=v/s)
      //=====================================================================================

      // Self-scaling (v=v/s)
      {
         test_ = "Self-scaling (v=v/s)";

         try {
            dres_   = lhs_ * rhs_;
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

         checkResults<MT>();
      }


      //=====================================================================================
      // Scaled multiplication (s*OP)
      //=====================================================================================

      // Scaled multiplication with the given matrix/vector
      {
         test_  = "Scaled multiplication with the given matrix/vector (s*OP)";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = scalar * ( lhs_ * rhs_ );
            sres_   = scalar * ( lhs_ * rhs_ );
            refres_ = scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = scalar * ( olhs_ * rhs_ );
            sres_   = scalar * ( olhs_ * rhs_ );
            refres_ = scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with evaluated matrix/vector
      {
         test_  = "Scaled multiplication with evaluated matrix/vector (s*OP)";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            sres_   = scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            refres_ = scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication (OP*s)
      //=====================================================================================

      // Scaled multiplication with the given matrix/vector
      {
         test_  = "Scaled multiplication with the given matrix/vector (OP*s)";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = ( lhs_ * rhs_ ) * scalar;
            sres_   = ( lhs_ * rhs_ ) * scalar;
            refres_ = ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = ( olhs_ * rhs_ ) * scalar;
            sres_   = ( olhs_ * rhs_ ) * scalar;
            refres_ = ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with evaluated matrix/vector
      {
         test_  = "Scaled multiplication with evaluated matrix/vector (OP*s)";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_ = ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            sres_ = ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication (OP/s)
      //=====================================================================================

      // Scaled multiplication with the given matrix/vector
      {
         test_  = "Scaled multiplication with the given matrix/vector (OP/s)";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = ( lhs_ * rhs_ ) / scalar;
            sres_   = ( lhs_ * rhs_ ) / scalar;
            refres_ = ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = ( olhs_ * rhs_ ) / scalar;
            sres_   = ( olhs_ * rhs_ ) / scalar;
            refres_ = ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with evaluated matrix/vector
      {
         test_  = "Scaled multiplication with evaluated matrix/vector (OP/s)";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_ = ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            sres_ = ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with addition assignment (s*OP)
      //=====================================================================================

      // Scaled multiplication with addition assignment with the given matrix/vector
      {
         test_  = "Scaled multiplication with addition assignment with the given matrix/vector (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * ( lhs_ * rhs_ );
            sres_   += scalar * ( lhs_ * rhs_ );
            refres_ += scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += scalar * ( olhs_ * rhs_ );
            sres_   += scalar * ( olhs_ * rhs_ );
            refres_ += scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with addition assignment with evaluated matrix/vector
      {
         test_  = "Scaled multiplication with addition assignment with evaluated matrix/vector (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   += scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            sres_   += scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with addition assignment (OP*s)
      //=====================================================================================

      // Scaled multiplication with addition assignment with the given matrix/vector
      {
         test_  = "Scaled multiplication with addition assignment with the given matrix/vector (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( lhs_ * rhs_ ) * scalar;
            sres_   += ( lhs_ * rhs_ ) * scalar;
            refres_ += ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += ( olhs_ * rhs_ ) * scalar;
            sres_   += ( olhs_ * rhs_ ) * scalar;
            refres_ += ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with addition assignment with evaluated matrix/vector
      {
         test_  = "Scaled multiplication with addition assignment with evaluated matrix/vector (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   += ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   += ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with addition assignment (OP/s)
      //=====================================================================================

      // Scaled multiplication with addition assignment with the given matrix/vector
      {
         test_  = "Scaled multiplication with addition assignment with the given matrix/vector (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( lhs_ * rhs_ ) / scalar;
            sres_   += ( lhs_ * rhs_ ) / scalar;
            refres_ += ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += ( olhs_ * rhs_ ) / scalar;
            sres_   += ( olhs_ * rhs_ ) / scalar;
            refres_ += ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with addition assignment with evaluated matrix/vector
      {
         test_  = "Scaled multiplication with addition assignment with evaluated matrix/vector (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   += ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   += ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled multiplication with subtraction assignment with the given matrix/vector
      {
         test_  = "Scaled multiplication with subtraction assignment with the given matrix/vector (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * ( lhs_ * rhs_ );
            sres_   -= scalar * ( lhs_ * rhs_ );
            refres_ -= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= scalar * ( olhs_ * rhs_ );
            sres_   -= scalar * ( olhs_ * rhs_ );
            refres_ -= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with subtraction assignment with evaluated matrix/vector
      {
         test_  = "Scaled multiplication with subtraction assignment with evaluated matrix/vector (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   -= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            sres_   -= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled multiplication with subtraction assignment with the given matrix/vector
      {
         test_  = "Scaled multiplication with subtraction assignment with the given matrix/vector (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( lhs_ * rhs_ ) * scalar;
            sres_   -= ( lhs_ * rhs_ ) * scalar;
            refres_ -= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= ( olhs_ * rhs_ ) * scalar;
            sres_   -= ( olhs_ * rhs_ ) * scalar;
            refres_ -= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with subtraction assignment with evaluated matrix/vector
      {
         test_  = "Scaled multiplication with subtraction assignment with evaluated matrix/vector (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   -= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled multiplication with subtraction assignment with the given matrix/vector
      {
         test_  = "Scaled multiplication with subtraction assignment with the given matrix/vector (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( lhs_ * rhs_ ) / scalar;
            sres_   -= ( lhs_ * rhs_ ) / scalar;
            refres_ -= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= ( olhs_ * rhs_ ) / scalar;
            sres_   -= ( olhs_ * rhs_ ) / scalar;
            refres_ -= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with subtraction assignment with evaluated matrix/vector
      {
         test_  = "Scaled multiplication with subtraction assignment with evaluated matrix/vector (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   -= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   -= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with multiplication assignment (s*OP)
      //=====================================================================================

      // Scaled multiplication with multiplication assignment with the given matrix/vector
      {
         test_  = "Scaled multiplication with multiplication assignment with the given matrix/vector (s*OP)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= scalar * ( lhs_ * rhs_ );
            sres_   *= scalar * ( lhs_ * rhs_ );
            refres_ *= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   *= scalar * ( olhs_ * rhs_ );
            sres_   *= scalar * ( olhs_ * rhs_ );
            refres_ *= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with multiplication assignment with evaluated matrix/vector
      {
         test_  = "Scaled multiplication with multiplication assignment with evaluated matrix/vector (s*OP)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   *= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ *= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   *= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            sres_   *= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            refres_ *= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with multiplication assignment (OP*s)
      //=====================================================================================

      // Scaled multiplication with multiplication assignment with the given matrix/vector
      {
         test_  = "Scaled multiplication with multiplication assignment with the given matrix/vector (OP*s)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= ( lhs_ * rhs_ ) * scalar;
            sres_   *= ( lhs_ * rhs_ ) * scalar;
            refres_ *= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   *= ( olhs_ * rhs_ ) * scalar;
            sres_   *= ( olhs_ * rhs_ ) * scalar;
            refres_ *= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with multiplication assignment with evaluated matrix/vector
      {
         test_  = "Scaled multiplication with multiplication assignment with evaluated matrix/vector (OP*s)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   *= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ *= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   *= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   *= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ *= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with multiplication assignment (OP/s)
      //=====================================================================================

      // Scaled multiplication with multiplication assignment with the given matrix/vector
      {
         test_  = "Scaled multiplication with multiplication assignment with the given matrix/vector (OP/s)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= ( lhs_ * rhs_ ) / scalar;
            sres_   *= ( lhs_ * rhs_ ) / scalar;
            refres_ *= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   *= ( olhs_ * rhs_ ) / scalar;
            sres_   *= ( olhs_ * rhs_ ) / scalar;
            refres_ *= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Scaled multiplication with multiplication assignment with evaluated matrix/vector
      {
         test_  = "Scaled multiplication with multiplication assignment with evaluated matrix/vector (OP/s)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   *= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ *= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   *= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   *= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ *= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Scaled multiplication with division assignment (s*OP)
      //=====================================================================================

      if( !blaze::IsUniform_v<MT> && blaze::isDivisor( lhs_ * rhs_ ) )
      {
         // Scaled multiplication with division assignment with the given matrix/vector
         {
            test_  = "Scaled multiplication with division assignment with the given matrix/vector (s*OP)";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= scalar * ( lhs_ * rhs_ );
               sres_   /= scalar * ( lhs_ * rhs_ );
               refres_ /= scalar * ( reflhs_ * refrhs_ );
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               dres_   /= scalar * ( olhs_ * rhs_ );
               sres_   /= scalar * ( olhs_ * rhs_ );
               refres_ /= scalar * ( reflhs_ * refrhs_ );
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }

         // Scaled multiplication with division assignment with evaluated matrix/vector
         {
            test_  = "Scaled multiplication with division assignment with evaluated matrix/vector (s*OP)";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
               sres_   /= scalar * ( eval( lhs_ ) * eval( rhs_ ) );
               refres_ /= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               dres_   /= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
               sres_   /= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
               refres_ /= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }
      }


      //=====================================================================================
      // Scaled multiplication with division assignment (OP*s)
      //=====================================================================================

      if( !blaze::IsUniform_v<MT> && blaze::isDivisor( lhs_ * rhs_ ) )
      {
         // Scaled multiplication with division assignment with the given matrix/vector
         {
            test_  = "Scaled multiplication with division assignment with the given matrix/vector (OP*s)";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= ( lhs_ * rhs_ ) * scalar;
               sres_   /= ( lhs_ * rhs_ ) * scalar;
               refres_ /= ( reflhs_ * refrhs_ ) * scalar;
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               dres_   /= ( olhs_ * rhs_ ) * scalar;
               sres_   /= ( olhs_ * rhs_ ) * scalar;
               refres_ /= ( reflhs_ * refrhs_ ) * scalar;
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }

         // Scaled multiplication with division assignment with evaluated matrix/vector
         {
            test_  = "Scaled multiplication with division assignment with evaluated matrix/vector (OP*s)";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
               sres_   /= ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
               refres_ /= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               dres_   /= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
               sres_   /= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
               refres_ /= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }
      }


      //=====================================================================================
      // Scaled multiplication with division assignment (OP/s)
      //=====================================================================================

      if( !blaze::IsUniform_v<MT> && blaze::isDivisor( ( lhs_ * rhs_ ) / scalar ) )
      {
         // Scaled multiplication with division assignment with the given matrix/vector
         {
            test_  = "Scaled multiplication with division assignment with the given matrix/vector (OP/s)";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= ( lhs_ * rhs_ ) / scalar;
               sres_   /= ( lhs_ * rhs_ ) / scalar;
               refres_ /= ( reflhs_ * refrhs_ ) / scalar;
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               dres_   /= ( olhs_ * rhs_ ) / scalar;
               sres_   /= ( olhs_ * rhs_ ) / scalar;
               refres_ /= ( reflhs_ * refrhs_ ) / scalar;
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }

         // Scaled multiplication with division assignment with evaluated matrix/vector
         {
            test_  = "Scaled multiplication with division assignment with evaluated matrix/vector (OP/s)";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
               sres_   /= ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
               refres_ /= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               dres_   /= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
               sres_   /= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
               refres_ /= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the transpose dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the transpose matrix/vector multiplication with plain assignment,
// addition assignment, subtraction assignment, multiplication assignment, and division
// assignment. In case any error resulting from the multiplication or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANS_OPERATION > 1 )
   {
      //=====================================================================================
      // Transpose multiplication
      //=====================================================================================

      // Transpose multiplication with the given matrix/vector
      {
         test_  = "Transpose multiplication with the given matrix/vector";
         error_ = "Failed multiplication operation";

         try {
            initTransposeResults();
            tdres_   = trans( lhs_ * rhs_ );
            tsres_   = trans( lhs_ * rhs_ );
            trefres_ = trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   = trans( olhs_ * rhs_ );
            tsres_   = trans( olhs_ * rhs_ );
            trefres_ = trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }

      // Transpose multiplication with evaluated matrix/vector
      {
         test_  = "Transpose multiplication with evaluated matrix/vector";
         error_ = "Failed multiplication operation";

         try {
            initTransposeResults();
            tdres_   = trans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_   = trans( eval( lhs_ ) * eval( rhs_ ) );
            trefres_ = trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   = trans( eval( olhs_ ) * eval( rhs_ ) );
            tsres_   = trans( eval( olhs_ ) * eval( rhs_ ) );
            trefres_ = trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }


      //=====================================================================================
      // Transpose multiplication with addition assignment
      //=====================================================================================

      // Transpose multiplication with addition assignment with the given matrix/vector
      {
         test_  = "Transpose multiplication with addition assignment with the given matrix/vector";
         error_ = "Failed addition assignment operation";

         try {
            initTransposeResults();
            tdres_   += trans( lhs_ * rhs_ );
            tsres_   += trans( lhs_ * rhs_ );
            trefres_ += trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   += trans( olhs_ * rhs_ );
            tsres_   += trans( olhs_ * rhs_ );
            trefres_ += trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }

      // Transpose multiplication with addition assignment with evaluated matrix/vector
      {
         test_  = "Transpose multiplication with addition assignment with evaluated matrix/vector";
         error_ = "Failed addition assignment operation";

         try {
            initTransposeResults();
            tdres_   += trans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_   += trans( eval( lhs_ ) * eval( rhs_ ) );
            trefres_ += trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   += trans( eval( olhs_ ) * eval( rhs_ ) );
            tsres_   += trans( eval( olhs_ ) * eval( rhs_ ) );
            trefres_ += trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }


      //=====================================================================================
      // Transpose multiplication with subtraction assignment
      //=====================================================================================

      // Transpose multiplication with subtraction assignment with the given matrix/vector
      {
         test_  = "Transpose multiplication with subtraction assignment with the given matrix/vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initTransposeResults();
            tdres_   -= trans( lhs_ * rhs_ );
            tsres_   -= trans( lhs_ * rhs_ );
            trefres_ -= trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   -= trans( olhs_ * rhs_ );
            tsres_   -= trans( olhs_ * rhs_ );
            trefres_ -= trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }

      // Transpose multiplication with subtraction assignment with evaluated matrix/vector
      {
         test_  = "Transpose multiplication with subtraction assignment with evaluated matrix/vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initTransposeResults();
            tdres_   -= trans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_   -= trans( eval( lhs_ ) * eval( rhs_ ) );
            trefres_ -= trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   -= trans( eval( olhs_ ) * eval( rhs_ ) );
            tsres_   -= trans( eval( olhs_ ) * eval( rhs_ ) );
            trefres_ -= trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }


      //=====================================================================================
      // Transpose multiplication with multiplication assignment
      //=====================================================================================

      // Transpose multiplication with multiplication assignment with the given matrix/vector
      {
         test_  = "Transpose multiplication with multiplication assignment with the given matrix/vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initTransposeResults();
            tdres_   *= trans( lhs_ * rhs_ );
            tsres_   *= trans( lhs_ * rhs_ );
            trefres_ *= trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   *= trans( olhs_ * rhs_ );
            tsres_   *= trans( olhs_ * rhs_ );
            trefres_ *= trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }

      // Transpose multiplication with multiplication assignment with evaluated matrix/vector
      {
         test_  = "Transpose multiplication with multiplication assignment with evaluated matrix/vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initTransposeResults();
            tdres_   *= trans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_   *= trans( eval( lhs_ ) * eval( rhs_ ) );
            trefres_ *= trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   *= trans( eval( olhs_ ) * eval( rhs_ ) );
            tsres_   *= trans( eval( olhs_ ) * eval( rhs_ ) );
            trefres_ *= trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }


      //=====================================================================================
      // Transpose multiplication with division assignment
      //=====================================================================================

      if( !blaze::IsUniform_v<MT> && blaze::isDivisor( lhs_ * rhs_ ) )
      {
         // Transpose multiplication with division assignment with the given matrix/vector
         {
            test_  = "Transpose multiplication with division assignment with the given matrix/vector";
            error_ = "Failed division assignment operation";

            try {
               initTransposeResults();
               tdres_   /= trans( lhs_ * rhs_ );
               tsres_   /= trans( lhs_ * rhs_ );
               trefres_ /= trans( reflhs_ * refrhs_ );
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkTransposeResults<MT>();

            try {
               initTransposeResults();
               tdres_   /= trans( olhs_ * rhs_ );
               tsres_   /= trans( olhs_ * rhs_ );
               trefres_ /= trans( reflhs_ * refrhs_ );
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkTransposeResults<TMT>();
         }

         // Transpose multiplication with division assignment with evaluated matrix/vector
         {
            test_  = "Transpose multiplication with division assignment with evaluated matrix/vector";
            error_ = "Failed division assignment operation";

            try {
               initTransposeResults();
               tdres_   /= trans( eval( lhs_ ) * eval( rhs_ ) );
               tsres_   /= trans( eval( lhs_ ) * eval( rhs_ ) );
               trefres_ /= trans( eval( reflhs_ ) * eval( refrhs_ ) );
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkTransposeResults<MT>();

            try {
               initTransposeResults();
               tdres_   /= trans( eval( olhs_ ) * eval( rhs_ ) );
               tsres_   /= trans( eval( olhs_ ) * eval( rhs_ ) );
               trefres_ /= trans( eval( reflhs_ ) * eval( refrhs_ ) );
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkTransposeResults<TMT>();
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the conjugate transpose dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the conjugate transpose matrix/vector multiplication with plain
// assignment, addition assignment, subtraction assignment, multiplication assignment,
// and division assignment. In case any error resulting from the multiplication or the
// subsequent assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testCTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION > 1 )
   {
      //=====================================================================================
      // Conjugate transpose multiplication
      //=====================================================================================

      // Conjugate transpose multiplication with the given matrix/vector
      {
         test_  = "Conjugate transpose multiplication with the given matrix/vector";
         error_ = "Failed multiplication operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( lhs_ * rhs_ );
            tsres_   = ctrans( lhs_ * rhs_ );
            trefres_ = ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   = ctrans( olhs_ * rhs_ );
            tsres_   = ctrans( olhs_ * rhs_ );
            trefres_ = ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }

      // Conjugate transpose multiplication with evaluated matrix/vector
      {
         test_  = "Conjugate transpose multiplication with evaluated matrix/vector";
         error_ = "Failed multiplication operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_   = ctrans( eval( lhs_ ) * eval( rhs_ ) );
            trefres_ = ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   = ctrans( eval( olhs_ ) * eval( rhs_ ) );
            tsres_   = ctrans( eval( olhs_ ) * eval( rhs_ ) );
            trefres_ = ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }


      //=====================================================================================
      // Conjugate transpose multiplication with addition assignment
      //=====================================================================================

      // Conjugate transpose multiplication with addition assignment with the given matrix/vector
      {
         test_  = "Conjugate transpose multiplication with addition assignment with the given matrix/vector";
         error_ = "Failed addition assignment operation";

         try {
            initTransposeResults();
            tdres_   += ctrans( lhs_ * rhs_ );
            tsres_   += ctrans( lhs_ * rhs_ );
            trefres_ += ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   += ctrans( olhs_ * rhs_ );
            tsres_   += ctrans( olhs_ * rhs_ );
            trefres_ += ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }

      // Conjugate transpose multiplication with addition assignment with evaluated matrix/vector
      {
         test_  = "Conjugate transpose multiplication with addition assignment with evaluated matrix/vector";
         error_ = "Failed addition assignment operation";

         try {
            initTransposeResults();
            tdres_   += ctrans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_   += ctrans( eval( lhs_ ) * eval( rhs_ ) );
            trefres_ += ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   += ctrans( eval( olhs_ ) * eval( rhs_ ) );
            tsres_   += ctrans( eval( olhs_ ) * eval( rhs_ ) );
            trefres_ += ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }


      //=====================================================================================
      // Conjugate transpose multiplication with subtraction assignment
      //=====================================================================================

      // Conjugate transpose multiplication with subtraction assignment with the given matrix/vector
      {
         test_  = "Conjugate transpose multiplication with subtraction assignment with the given matrix/vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initTransposeResults();
            tdres_   -= ctrans( lhs_ * rhs_ );
            tsres_   -= ctrans( lhs_ * rhs_ );
            trefres_ -= ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   -= ctrans( olhs_ * rhs_ );
            tsres_   -= ctrans( olhs_ * rhs_ );
            trefres_ -= ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }

      // Conjugate transpose multiplication with subtraction assignment with evaluated matrix/vector
      {
         test_  = "Conjugate transpose multiplication with subtraction assignment with evaluated matrix/vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initTransposeResults();
            tdres_   -= ctrans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_   -= ctrans( eval( lhs_ ) * eval( rhs_ ) );
            trefres_ -= ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   -= ctrans( eval( olhs_ ) * eval( rhs_ ) );
            tsres_   -= ctrans( eval( olhs_ ) * eval( rhs_ ) );
            trefres_ -= ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }


      //=====================================================================================
      // Conjugate transpose multiplication with multiplication assignment
      //=====================================================================================

      // Conjugate transpose multiplication with multiplication assignment with the given matrix/vector
      {
         test_  = "Conjugate transpose multiplication with multiplication assignment with the given matrix/vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initTransposeResults();
            tdres_   *= ctrans( lhs_ * rhs_ );
            tsres_   *= ctrans( lhs_ * rhs_ );
            trefres_ *= ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   *= ctrans( olhs_ * rhs_ );
            tsres_   *= ctrans( olhs_ * rhs_ );
            trefres_ *= ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }

      // Conjugate transpose multiplication with multiplication assignment with evaluated matrix/vector
      {
         test_  = "Conjugate transpose multiplication with multiplication assignment with evaluated matrix/vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initTransposeResults();
            tdres_   *= ctrans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_   *= ctrans( eval( lhs_ ) * eval( rhs_ ) );
            trefres_ *= ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_   *= ctrans( eval( olhs_ ) * eval( rhs_ ) );
            tsres_   *= ctrans( eval( olhs_ ) * eval( rhs_ ) );
            trefres_ *= ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkTransposeResults<TMT>();
      }


      //=====================================================================================
      // Conjugate transpose multiplication with division assignment
      //=====================================================================================

      if( !blaze::IsUniform_v<MT> && blaze::isDivisor( lhs_ * rhs_ ) )
      {
         // Conjugate transpose multiplication with division assignment with the given matrix/vector
         {
            test_  = "Conjugate transpose multiplication with division assignment with the given matrix/vector";
            error_ = "Failed division assignment operation";

            try {
               initTransposeResults();
               tdres_   /= ctrans( lhs_ * rhs_ );
               tsres_   /= ctrans( lhs_ * rhs_ );
               trefres_ /= ctrans( reflhs_ * refrhs_ );
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkTransposeResults<MT>();

            try {
               initTransposeResults();
               tdres_   /= ctrans( olhs_ * rhs_ );
               tsres_   /= ctrans( olhs_ * rhs_ );
               trefres_ /= ctrans( reflhs_ * refrhs_ );
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkTransposeResults<TMT>();
         }

         // Conjugate transpose multiplication with division assignment with evaluated matrix/vector
         {
            test_  = "Conjugate transpose multiplication with division assignment with evaluated matrix/vector";
            error_ = "Failed division assignment operation";

            try {
               initTransposeResults();
               tdres_   /= ctrans( eval( lhs_ ) * eval( rhs_ ) );
               tsres_   /= ctrans( eval( lhs_ ) * eval( rhs_ ) );
               trefres_ /= ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkTransposeResults<MT>();

            try {
               initTransposeResults();
               tdres_   /= ctrans( eval( olhs_ ) * eval( rhs_ ) );
               tsres_   /= ctrans( eval( olhs_ ) * eval( rhs_ ) );
               trefres_ /= ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkTransposeResults<TMT>();
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the abs dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the abs matrix/vector multiplication with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the multiplication or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testAbsOperation()
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
/*!\brief Testing the conjugate dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the conjugate matrix/vector multiplication with plain assignment,
// addition assignment, subtraction assignment, multiplication assignment, and division
// assignment. In case any error resulting from the multiplication or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testConjOperation()
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
/*!\brief Testing the \a real dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the \a real matrix/vector multiplication with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the multiplication or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testRealOperation()
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
/*!\brief Testing the \a imag dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the \a imag matrix/vector multiplication with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the multiplication or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testImagOperation()
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
/*!\brief Testing the evaluated dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the evaluated matrix/vector multiplication with plain assignment,
// addition assignment, subtraction assignment, multiplication assignment, and division
// assignment. In case any error resulting from the multiplication or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testEvalOperation()
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
/*!\brief Testing the serialized dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the serialized matrix/vector multiplication with plain assignment,
// addition assignment, subtraction assignment, multiplication assignment, and division
// assignment. In case any error resulting from the multiplication or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testSerialOperation()
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
/*!\brief Testing the non-aliased dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the non-aliased matrix/vector multiplication with plain assignment,
// addition assignment, subtraction assignment, multiplication assignment, and division
// assignment. In case any error resulting from the multiplication or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testNoAliasOperation()
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
/*!\brief Testing the non-SIMD dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the non-SIMD matrix/vector multiplication with plain assignment,
// addition assignment, subtraction assignment, multiplication assignment, and division
// assignment. In case any error resulting from the multiplication or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testNoSIMDOperation()
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
/*!\brief Testing the subvector-wise dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the subvector-wise matrix/vector multiplication with plain assignment,
// addition assignment, subtraction assignment, multiplication assignment, and division
// assignment. In case any error resulting from the multiplication or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testSubvectorOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION > 1 )
   {
      if( lhs_.rows() == 0UL )
         return;


      //=====================================================================================
      // Subvector-wise multiplication
      //=====================================================================================

      // Subvector-wise multiplication with the given matrix/vector
      {
         test_  = "Subvector-wise multiplication with the given matrix/vector";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.rows() - index );
               subvector( dres_  , index, size ) = subvector( lhs_ * rhs_      , index, size );
               subvector( sres_  , index, size ) = subvector( lhs_ * rhs_      , index, size );
               subvector( refres_, index, size ) = subvector( reflhs_ * refrhs_, index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<olhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, olhs_.rows() - index );
               subvector( dres_  , index, size ) = subvector( olhs_ * rhs_     , index, size );
               subvector( sres_  , index, size ) = subvector( olhs_ * rhs_     , index, size );
               subvector( refres_, index, size ) = subvector( reflhs_ * refrhs_, index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Subvector-wise multiplication with evaluated matrix/vector
      {
         test_  = "Subvector-wise multiplication with evaluated matrix/vector";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.rows() - index );
               subvector( dres_  , index, size ) = subvector( eval( lhs_ ) * eval( rhs_ )      , index, size );
               subvector( sres_  , index, size ) = subvector( eval( lhs_ ) * eval( rhs_ )      , index, size );
               subvector( refres_, index, size ) = subvector( eval( reflhs_ ) * eval( refrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<olhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, olhs_.rows() - index );
               subvector( dres_  , index, size ) = subvector( eval( olhs_ ) * eval( rhs_ )     , index, size );
               subvector( sres_  , index, size ) = subvector( eval( olhs_ ) * eval( rhs_ )     , index, size );
               subvector( refres_, index, size ) = subvector( eval( reflhs_ ) * eval( refrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Subvector-wise multiplication with addition assignment
      //=====================================================================================

      // Subvector-wise multiplication with addition assignment with the given matrix/vector
      {
         test_  = "Subvector-wise multiplication with addition assignment the given matrix/vector";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.rows() - index );
               subvector( dres_  , index, size ) += subvector( lhs_ * rhs_      , index, size );
               subvector( sres_  , index, size ) += subvector( lhs_ * rhs_      , index, size );
               subvector( refres_, index, size ) += subvector( reflhs_ * refrhs_, index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<olhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, olhs_.rows() - index );
               subvector( dres_  , index, size ) += subvector( olhs_ * rhs_     , index, size );
               subvector( sres_  , index, size ) += subvector( olhs_ * rhs_     , index, size );
               subvector( refres_, index, size ) += subvector( reflhs_ * refrhs_, index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Subvector-wise multiplication with addition assignment with evaluated matrix/vector
      {
         test_  = "Subvector-wise multiplication with addition assignment with evaluated matrix/vector";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.rows() - index );
               subvector( dres_  , index, size ) += subvector( eval( lhs_ ) * eval( rhs_ )      , index, size );
               subvector( sres_  , index, size ) += subvector( eval( lhs_ ) * eval( rhs_ )      , index, size );
               subvector( refres_, index, size ) += subvector( eval( reflhs_ ) * eval( refrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<olhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, olhs_.rows() - index );
               subvector( dres_  , index, size ) += subvector( eval( olhs_ ) * eval( rhs_ )     , index, size );
               subvector( sres_  , index, size ) += subvector( eval( olhs_ ) * eval( rhs_ )     , index, size );
               subvector( refres_, index, size ) += subvector( eval( reflhs_ ) * eval( refrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Subvector-wise multiplication with subtraction assignment
      //=====================================================================================

      // Subvector-wise multiplication with addition assignment with the given matrix/vector
      {
         test_  = "Subvector-wise multiplication with subtraction assignment the given matrix/vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.rows() - index );
               subvector( dres_  , index, size ) -= subvector( lhs_ * rhs_      , index, size );
               subvector( sres_  , index, size ) -= subvector( lhs_ * rhs_      , index, size );
               subvector( refres_, index, size ) -= subvector( reflhs_ * refrhs_, index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<olhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, olhs_.rows() - index );
               subvector( dres_  , index, size ) -= subvector( olhs_ * rhs_     , index, size );
               subvector( sres_  , index, size ) -= subvector( olhs_ * rhs_     , index, size );
               subvector( refres_, index, size ) -= subvector( reflhs_ * refrhs_, index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Subvector-wise multiplication with subtraction assignment with evaluated matrix/vector
      {
         test_  = "Subvector-wise multiplication with subtraction assignment with evaluated matrix/vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.rows() - index );
               subvector( dres_  , index, size ) -= subvector( eval( lhs_ ) * eval( rhs_ )      , index, size );
               subvector( sres_  , index, size ) -= subvector( eval( lhs_ ) * eval( rhs_ )      , index, size );
               subvector( refres_, index, size ) -= subvector( eval( reflhs_ ) * eval( refrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<olhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, olhs_.rows() - index );
               subvector( dres_  , index, size ) -= subvector( eval( olhs_ ) * eval( rhs_ )     , index, size );
               subvector( sres_  , index, size ) -= subvector( eval( olhs_ ) * eval( rhs_ )     , index, size );
               subvector( refres_, index, size ) -= subvector( eval( reflhs_ ) * eval( refrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Subvector-wise multiplication with multiplication assignment
      //=====================================================================================

      // Subvector-wise multiplication with multiplication assignment with the given matrix/vector
      {
         test_  = "Subvector-wise multiplication with multiplication assignment the given matrix/vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.rows() - index );
               subvector( dres_  , index, size ) *= subvector( lhs_ * rhs_      , index, size );
               subvector( sres_  , index, size ) *= subvector( lhs_ * rhs_      , index, size );
               subvector( refres_, index, size ) *= subvector( reflhs_ * refrhs_, index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<olhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, olhs_.rows() - index );
               subvector( dres_  , index, size ) *= subvector( olhs_ * rhs_     , index, size );
               subvector( sres_  , index, size ) *= subvector( olhs_ * rhs_     , index, size );
               subvector( refres_, index, size ) *= subvector( reflhs_ * refrhs_, index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Subvector-wise multiplication with multiplication assignment with evaluated matrix/vector
      {
         test_  = "Subvector-wise multiplication with multiplication assignment with evaluated matrix/vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.rows() - index );
               subvector( dres_  , index, size ) *= subvector( eval( lhs_ ) * eval( rhs_ )      , index, size );
               subvector( sres_  , index, size ) *= subvector( eval( lhs_ ) * eval( rhs_ )      , index, size );
               subvector( refres_, index, size ) *= subvector( eval( reflhs_ ) * eval( refrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<olhs_.rows(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, olhs_.rows() - index );
               subvector( dres_  , index, size ) *= subvector( eval( olhs_ ) * eval( rhs_ )     , index, size );
               subvector( sres_  , index, size ) *= subvector( eval( olhs_ ) * eval( rhs_ )     , index, size );
               subvector( refres_, index, size ) *= subvector( eval( reflhs_ ) * eval( refrhs_ ), index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Subvector-wise multiplication with division assignment
      //=====================================================================================

      if( !blaze::IsUniform_v<MT> )
      {
         // Subvector-wise multiplication with division assignment with the given matrix/vector
         {
            test_  = "Subvector-wise multiplication with division assignment the given matrix/vector";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               for( size_t index=0UL, size=0UL; index<lhs_.rows(); index+=size ) {
                  size = blaze::rand<size_t>( 1UL, lhs_.rows() - index );
                  if( !blaze::isDivisor( subvector( lhs_ * rhs_, index, size ) ) ) continue;
                  subvector( dres_  , index, size ) /= subvector( lhs_ * rhs_      , index, size );
                  subvector( sres_  , index, size ) /= subvector( lhs_ * rhs_      , index, size );
                  subvector( refres_, index, size ) /= subvector( reflhs_ * refrhs_, index, size );
               }
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               for( size_t index=0UL, size=0UL; index<olhs_.rows(); index+=size ) {
                  size = blaze::rand<size_t>( 1UL, olhs_.rows() - index );
                  if( !blaze::isDivisor( subvector( olhs_ * rhs_, index, size ) ) ) continue;
                  subvector( dres_  , index, size ) /= subvector( olhs_ * rhs_     , index, size );
                  subvector( sres_  , index, size ) /= subvector( olhs_ * rhs_     , index, size );
                  subvector( refres_, index, size ) /= subvector( reflhs_ * refrhs_, index, size );
               }
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }

         // Subvector-wise multiplication with division assignment with evaluated matrix/vector
         {
            test_  = "Subvector-wise multiplication with division assignment with evaluated matrix/vector";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               for( size_t index=0UL, size=0UL; index<lhs_.rows(); index+=size ) {
                  size = blaze::rand<size_t>( 1UL, lhs_.rows() - index );
                  if( !blaze::isDivisor( subvector( lhs_ * rhs_, index, size ) ) ) continue;
                  subvector( dres_  , index, size ) /= subvector( eval( lhs_ ) * eval( rhs_ )      , index, size );
                  subvector( sres_  , index, size ) /= subvector( eval( lhs_ ) * eval( rhs_ )      , index, size );
                  subvector( refres_, index, size ) /= subvector( eval( reflhs_ ) * eval( refrhs_ ), index, size );
               }
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               for( size_t index=0UL, size=0UL; index<olhs_.rows(); index+=size ) {
                  size = blaze::rand<size_t>( 1UL, olhs_.rows() - index );
                  if( !blaze::isDivisor( subvector( olhs_ * rhs_, index, size ) ) ) continue;
                  subvector( dres_  , index, size ) /= subvector( eval( olhs_ ) * eval( rhs_ )     , index, size );
                  subvector( sres_  , index, size ) /= subvector( eval( olhs_ ) * eval( rhs_ )     , index, size );
                  subvector( refres_, index, size ) /= subvector( eval( reflhs_ ) * eval( refrhs_ ), index, size );
               }
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      try {
         auto sv = subvector( lhs_ * rhs_, 1UL, lhs_.rows() );

         std::ostringstream oss;
         oss << " Test: Subvector construction\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense matrix type:\n"
             << "     " << typeid( MT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
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
         auto sv = subvector( lhs_ * rhs_, lhs_.rows(), 1UL );

         std::ostringstream oss;
         oss << " Test: Subvector construction\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense matrix type:\n"
             << "     " << typeid( MT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
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
/*!\brief Skipping the subvector-wise dense matrix/sparse vector multiplication.
//
// \return void
//
// This function is called in case the subvector-wise matrix/vector multiplication operation is
// not available for the given types \a MT and \a VT.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testSubvectorOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the elements-wise dense matrix/sparse vector multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the elements-wise matrix/vector multiplication with plain assignment,
// addition assignment, subtraction assignment, multiplication assignment, and division
// assignment. In case any error resulting from the multiplication or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testElementsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_ELEMENTS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ELEMENTS_OPERATION > 1 )
   {
      if( lhs_.rows() == 0UL )
         return;


      std::vector<size_t> indices( lhs_.rows() );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );


      //=====================================================================================
      // Elements-wise multiplication
      //=====================================================================================

      // Elements-wise multiplication with the given matrix/vector
      {
         test_  = "Elements-wise multiplication with the given matrix/vector";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) = elements( lhs_ * rhs_      , &indices[index], n );
               elements( sres_  , &indices[index], n ) = elements( lhs_ * rhs_      , &indices[index], n );
               elements( refres_, &indices[index], n ) = elements( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) = elements( olhs_ * rhs_     , &indices[index], n );
               elements( sres_  , &indices[index], n ) = elements( olhs_ * rhs_     , &indices[index], n );
               elements( refres_, &indices[index], n ) = elements( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Elements-wise multiplication with evaluated matrix/vector
      {
         test_  = "Elements-wise multiplication with evaluated matrix/vector";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) = elements( eval( lhs_ ) * eval( rhs_ )      , &indices[index], n );
               elements( sres_  , &indices[index], n ) = elements( eval( lhs_ ) * eval( rhs_ )      , &indices[index], n );
               elements( refres_, &indices[index], n ) = elements( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) = elements( eval( olhs_ ) * eval( rhs_ )     , &indices[index], n );
               elements( sres_  , &indices[index], n ) = elements( eval( olhs_ ) * eval( rhs_ )     , &indices[index], n );
               elements( refres_, &indices[index], n ) = elements( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Elements-wise multiplication with addition assignment
      //=====================================================================================

      // Elements-wise multiplication with addition assignment with the given matrix/vector
      {
         test_  = "Elements-wise multiplication with addition assignment the given matrix/vector";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) += elements( lhs_ * rhs_      , &indices[index], n );
               elements( sres_  , &indices[index], n ) += elements( lhs_ * rhs_      , &indices[index], n );
               elements( refres_, &indices[index], n ) += elements( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) += elements( olhs_ * rhs_     , &indices[index], n );
               elements( sres_  , &indices[index], n ) += elements( olhs_ * rhs_     , &indices[index], n );
               elements( refres_, &indices[index], n ) += elements( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Elements-wise multiplication with addition assignment with evaluated matrix/vector
      {
         test_  = "Elements-wise multiplication with addition assignment with evaluated matrix/vector";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) += elements( eval( lhs_ ) * eval( rhs_ )      , &indices[index], n );
               elements( sres_  , &indices[index], n ) += elements( eval( lhs_ ) * eval( rhs_ )      , &indices[index], n );
               elements( refres_, &indices[index], n ) += elements( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) += elements( eval( olhs_ ) * eval( rhs_ )     , &indices[index], n );
               elements( sres_  , &indices[index], n ) += elements( eval( olhs_ ) * eval( rhs_ )     , &indices[index], n );
               elements( refres_, &indices[index], n ) += elements( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Elements-wise multiplication with subtraction assignment
      //=====================================================================================

      // Elements-wise multiplication with subtraction assignment with the given matrix/vector
      {
         test_  = "Elements-wise multiplication with subtraction assignment the given matrix/vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) -= elements( lhs_ * rhs_      , &indices[index], n );
               elements( sres_  , &indices[index], n ) -= elements( lhs_ * rhs_      , &indices[index], n );
               elements( refres_, &indices[index], n ) -= elements( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) -= elements( olhs_ * rhs_     , &indices[index], n );
               elements( sres_  , &indices[index], n ) -= elements( olhs_ * rhs_     , &indices[index], n );
               elements( refres_, &indices[index], n ) -= elements( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Elements-wise multiplication with subtraction assignment with evaluated matrix/vector
      {
         test_  = "Elements-wise multiplication with subtraction assignment with evaluated matrix/vector";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) -= elements( eval( lhs_ ) * eval( rhs_ )      , &indices[index], n );
               elements( sres_  , &indices[index], n ) -= elements( eval( lhs_ ) * eval( rhs_ )      , &indices[index], n );
               elements( refres_, &indices[index], n ) -= elements( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) -= elements( eval( olhs_ ) * eval( rhs_ )     , &indices[index], n );
               elements( sres_  , &indices[index], n ) -= elements( eval( olhs_ ) * eval( rhs_ )     , &indices[index], n );
               elements( refres_, &indices[index], n ) -= elements( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Elements-wise multiplication with multiplication assignment
      //=====================================================================================

      // Elements-wise multiplication with multiplication assignment with the given matrix/vector
      {
         test_  = "Elements-wise multiplication with multiplication assignment the given matrix/vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) *= elements( lhs_ * rhs_      , &indices[index], n );
               elements( sres_  , &indices[index], n ) *= elements( lhs_ * rhs_      , &indices[index], n );
               elements( refres_, &indices[index], n ) *= elements( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) *= elements( olhs_ * rhs_     , &indices[index], n );
               elements( sres_  , &indices[index], n ) *= elements( olhs_ * rhs_     , &indices[index], n );
               elements( refres_, &indices[index], n ) *= elements( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Elements-wise multiplication with multiplication assignment with evaluated matrix/vector
      {
         test_  = "Elements-wise multiplication with multiplication assignment with evaluated matrix/vector";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) *= elements( eval( lhs_ ) * eval( rhs_ )      , &indices[index], n );
               elements( sres_  , &indices[index], n ) *= elements( eval( lhs_ ) * eval( rhs_ )      , &indices[index], n );
               elements( refres_, &indices[index], n ) *= elements( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) *= elements( eval( olhs_ ) * eval( rhs_ )     , &indices[index], n );
               elements( sres_  , &indices[index], n ) *= elements( eval( olhs_ ) * eval( rhs_ )     , &indices[index], n );
               elements( refres_, &indices[index], n ) *= elements( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }


      //=====================================================================================
      // Elements-wise multiplication with division assignment
      //=====================================================================================

      if( !blaze::IsUniform_v<MT> )
      {
         // Elements-wise multiplication with division assignment with the given matrix/vector
         {
            test_  = "Elements-wise multiplication with division assignment the given matrix/vector";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
                  n = blaze::rand<size_t>( 1UL, indices.size() - index );
                  if( !blaze::isDivisor( elements( lhs_ * rhs_, &indices[index], n ) ) ) continue;
                  elements( dres_  , &indices[index], n ) /= elements( lhs_ * rhs_      , &indices[index], n );
                  elements( sres_  , &indices[index], n ) /= elements( lhs_ * rhs_      , &indices[index], n );
                  elements( refres_, &indices[index], n ) /= elements( reflhs_ * refrhs_, &indices[index], n );
               }
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
                  n = blaze::rand<size_t>( 1UL, indices.size() - index );
                  if( !blaze::isDivisor( elements( olhs_ * rhs_, &indices[index], n ) ) ) continue;
                  elements( dres_  , &indices[index], n ) /= elements( olhs_ * rhs_     , &indices[index], n );
                  elements( sres_  , &indices[index], n ) /= elements( olhs_ * rhs_     , &indices[index], n );
                  elements( refres_, &indices[index], n ) /= elements( reflhs_ * refrhs_, &indices[index], n );
               }
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }

         // Elements-wise multiplication with division assignment with evaluated matrix/vector
         {
            test_  = "Elements-wise multiplication with division assignment with evaluated matrix/vector";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
                  n = blaze::rand<size_t>( 1UL, indices.size() - index );
                  if( !blaze::isDivisor( elements( lhs_ * rhs_, &indices[index], n ) ) ) continue;
                  elements( dres_  , &indices[index], n ) /= elements( eval( lhs_ ) * eval( rhs_ )      , &indices[index], n );
                  elements( sres_  , &indices[index], n ) /= elements( eval( lhs_ ) * eval( rhs_ )      , &indices[index], n );
                  elements( refres_, &indices[index], n ) /= elements( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
               }
            }
            catch( std::exception& ex ) {
               convertException<MT>( ex );
            }

            checkResults<MT>();

            try {
               initResults();
               for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
                  n = blaze::rand<size_t>( 1UL, indices.size() - index );
                  if( !blaze::isDivisor( elements( olhs_ * rhs_, &indices[index], n ) ) ) continue;
                  elements( dres_  , &indices[index], n ) /= elements( eval( olhs_ ) * eval( rhs_ )     , &indices[index], n );
                  elements( sres_  , &indices[index], n ) /= elements( eval( olhs_ ) * eval( rhs_ )     , &indices[index], n );
                  elements( refres_, &indices[index], n ) /= elements( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
               }
            }
            catch( std::exception& ex ) {
               convertException<TMT>( ex );
            }

            checkResults<TMT>();
         }
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      try {
         auto e = elements( lhs_ * rhs_, blaze::index_sequence<128>() );

         std::ostringstream oss;
         oss << " Test: Elements construction (index_sequence)\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense matrix type:\n"
             << "     " << typeid( MT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
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
         auto e = elements( lhs_ * rhs_, { lhs_.rows() } );

         std::ostringstream oss;
         oss << " Test: Elements construction (initializer_list)\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense matrix type:\n"
             << "     " << typeid( MT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
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
         auto e = elements( lhs_ * rhs_, [index=lhs_.rows()]( size_t ){ return index; }, 1UL );

         std::ostringstream oss;
         oss << " Test: Elements construction (lambda)\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense matrix type:\n"
             << "     " << typeid( MT ).name() << "\n"
             << "   Right-hand side sparse vector type:\n"
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
/*!\brief Skipping the elements-wise dense matrix/sparse vector multiplication.
//
// \return void
//
// This function is called in case the elements-wise matrix/vector multiplication operation is
// not available for the given types \a MT and \a VT.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::testElementsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the customized dense matrix/sparse vector multiplication.
//
// \param op The custom operation to be tested.
// \param name The human-readable name of the operation.
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the matrix/vector multiplication with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment
// in combination with a custom operation. In case any error resulting from the multiplication
// or the subsequent assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
template< typename OP >  // Type of the custom operation
void OperationTest<MT,VT>::testCustomOperation( OP op, const std::string& name )
{
   //=====================================================================================
   // Customized multiplication
   //=====================================================================================

   // Customized multiplication with the given matrix/vector
   {
      test_  = "Customized multiplication with the given matrix/vector";
      error_ = "Failed multiplication operation";

      try {
         initResults();
         dres_   = op( lhs_ * rhs_ );
         sres_   = op( lhs_ * rhs_ );
         refres_ = op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   = op( olhs_ * rhs_ );
         sres_   = op( olhs_ * rhs_ );
         refres_ = op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<TMT>( ex );
      }

      checkResults<TMT>();
   }

   // Customized multiplication with evaluated matrix/vector
   {
      test_  = "Customized multiplication with evaluated matrix/vector (" + name + ")";
      error_ = "Failed multiplication operation";

      try {
         initResults();
         dres_   = op( eval( lhs_ ) * eval( rhs_ ) );
         sres_   = op( eval( lhs_ ) * eval( rhs_ ) );
         refres_ = op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   = op( eval( olhs_ ) * eval( rhs_ ) );
         sres_   = op( eval( olhs_ ) * eval( rhs_ ) );
         refres_ = op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TMT>( ex );
      }

      checkResults<TMT>();
   }


   //=====================================================================================
   // Customized multiplication with addition assignment
   //=====================================================================================

   // Customized multiplication with addition assignment with the given matrix/vector
   {
      test_  = "Customized multiplication with addition assignment with the given matrix/vector (" + name + ")";
      error_ = "Failed addition assignment operation";

      try {
         initResults();
         dres_   += op( lhs_ * rhs_ );
         sres_   += op( lhs_ * rhs_ );
         refres_ += op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   += op( olhs_ * rhs_ );
         sres_   += op( olhs_ * rhs_ );
         refres_ += op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<TMT>( ex );
      }

      checkResults<TMT>();
   }

   // Customized multiplication with addition assignment with evaluated matrix/vector
   {
      test_  = "Customized multiplication with addition assignment with evaluated matrix/vector (" + name + ")";
      error_ = "Failed addition assignment operation";

      try {
         initResults();
         dres_   += op( eval( lhs_ ) * eval( rhs_ ) );
         sres_   += op( eval( lhs_ ) * eval( rhs_ ) );
         refres_ += op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   += op( eval( olhs_ ) * eval( rhs_ ) );
         sres_   += op( eval( olhs_ ) * eval( rhs_ ) );
         refres_ += op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TMT>( ex );
      }

      checkResults<TMT>();
   }


   //=====================================================================================
   // Customized multiplication with subtraction assignment
   //=====================================================================================

   // Customized multiplication with subtraction assignment with the given matrix/vector
   {
      test_  = "Customized multiplication with subtraction assignment with the given matrix/vector (" + name + ")";
      error_ = "Failed subtraction assignment operation";

      try {
         initResults();
         dres_   -= op( lhs_ * rhs_ );
         sres_   -= op( lhs_ * rhs_ );
         refres_ -= op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   -= op( olhs_ * rhs_ );
         sres_   -= op( olhs_ * rhs_ );
         refres_ -= op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<TMT>( ex );
      }

      checkResults<TMT>();
   }

   // Customized multiplication with subtraction assignment with evaluated matrix/vector
   {
      test_  = "Customized multiplication with subtraction assignment with evaluated matrix/vector (" + name + ")";
      error_ = "Failed subtraction assignment operation";

      try {
         initResults();
         dres_   -= op( eval( lhs_ ) * eval( rhs_ ) );
         sres_   -= op( eval( lhs_ ) * eval( rhs_ ) );
         refres_ -= op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   -= op( eval( olhs_ ) * eval( rhs_ ) );
         sres_   -= op( eval( olhs_ ) * eval( rhs_ ) );
         refres_ -= op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TMT>( ex );
      }

      checkResults<TMT>();
   }


   //=====================================================================================
   // Customized multiplication with multiplication assignment
   //=====================================================================================

   // Customized multiplication with multiplication assignment with the given matrix/vector
   {
      test_  = "Customized multiplication with multiplication assignment with the given matrix/vector (" + name + ")";
      error_ = "Failed multiplication assignment operation";

      try {
         initResults();
         dres_   *= op( lhs_ * rhs_ );
         sres_   *= op( lhs_ * rhs_ );
         refres_ *= op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   *= op( olhs_ * rhs_ );
         sres_   *= op( olhs_ * rhs_ );
         refres_ *= op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<TMT>( ex );
      }

      checkResults<TMT>();
   }

   // Customized multiplication with multiplication assignment with evaluated matrix/vector
   {
      test_  = "Customized multiplication with multiplication assignment with evaluated matrix/vector (" + name + ")";
      error_ = "Failed multiplication assignment operation";

      try {
         initResults();
         dres_   *= op( eval( lhs_ ) * eval( rhs_ ) );
         sres_   *= op( eval( lhs_ ) * eval( rhs_ ) );
         refres_ *= op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   *= op( eval( olhs_ ) * eval( rhs_ ) );
         sres_   *= op( eval( olhs_ ) * eval( rhs_ ) );
         refres_ *= op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TMT>( ex );
      }

      checkResults<TMT>();
   }


   //=====================================================================================
   // Customized multiplication with division assignment
   //=====================================================================================

   if( !blaze::IsUniform_v<VT> && !blaze::IsUniform_v<MT> && blaze::isDivisor( op( lhs_ * rhs_ ) ) )
   {
      // Customized multiplication with division assignment with the given matrix/vector
      {
         test_  = "Customized multiplication with division assignment with the given matrix/vector (" + name + ")";
         error_ = "Failed division assignment operation";

         try {
            initResults();
            dres_   /= op( lhs_ * rhs_ );
            sres_   /= op( lhs_ * rhs_ );
            refres_ /= op( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   /= op( olhs_ * rhs_ );
            sres_   /= op( olhs_ * rhs_ );
            refres_ /= op( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }

      // Customized multiplication with division assignment with evaluated matrix/vector
      {
         test_  = "Customized multiplication with division assignment with evaluated matrix/vector (" + name + ")";
         error_ = "Failed division assignment operation";

         try {
            initResults();
            dres_   /= op( eval( lhs_ ) * eval( rhs_ ) );
            sres_   /= op( eval( lhs_ ) * eval( rhs_ ) );
            refres_ /= op( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   /= op( eval( olhs_ ) * eval( rhs_ ) );
            sres_   /= op( eval( olhs_ ) * eval( rhs_ ) );
            refres_ /= op( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TMT>( ex );
         }

         checkResults<TMT>();
      }
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
// The template argument \a LT indicates the types of the left-hand side operand used for
// the computations.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
template< typename LT >  // Type of the left-hand side operand
void OperationTest<MT,VT>::checkResults()
{
   using blaze::IsRowMajorMatrix;

   if( !isEqual( dres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Result:\n" << dres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( sres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
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
// results. The template argument \a LT indicates the types of the left-hand side operand
// used for the computations.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
template< typename LT >  // Type of the left-hand side operand
void OperationTest<MT,VT>::checkTransposeResults()
{
   using blaze::IsRowMajorMatrix;

   if( !isEqual( tdres_, trefres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Transpose result:\n" << tdres_ << "\n"
          << "   Expected transpose result:\n" << trefres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( tsres_, trefres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side sparse vector type:\n"
          << "     " << typeid( VT ).name() << "\n"
          << "   Transpose result:\n" << tsres_ << "\n"
          << "   Expected transpose result:\n" << trefres_ << "\n";
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
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::initResults()
{
   const blaze::UnderlyingBuiltin_t<DRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<DRE> max( randmax );

   resize( dres_, rows( lhs_ ) );
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
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void OperationTest<MT,VT>::initTransposeResults()
{
   const blaze::UnderlyingBuiltin_t<TDRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<TDRE> max( randmax );

   resize( tdres_, rows( lhs_ ) );
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
// test. The template argument \a LT indicates the types of the left-hand side operand used for
// used for the computations.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
template< typename LT >  // Type of the left-hand side operand
void OperationTest<MT,VT>::convertException( const std::exception& ex )
{
   using blaze::IsRowMajorMatrix;

   std::ostringstream oss;
   oss << " Test : " << test_ << "\n"
       << " Error: " << error_ << "\n"
       << " Details:\n"
       << "   Random seed = " << blaze::getSeed() << "\n"
       << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
       << "     " << typeid( LT ).name() << "\n"
       << "   Right-hand side sparse vector type:\n"
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
/*!\brief Testing the matrix/vector multiplication between two specific types.
//
// \param creator1 The creator for the left-hand side matrix.
// \param creator2 The creator for the right-hand side vector.
// \return void
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
void runTest( const Creator<MT>& creator1, const Creator<VT>& creator2 )
{
#if BLAZETEST_MATHTEST_TEST_MULTIPLICATION
   if( BLAZETEST_MATHTEST_TEST_MULTIPLICATION > 1 )
   {
      for( size_t rep=0UL; rep<BLAZETEST_REPETITIONS; ++rep ) {
         OperationTest<MT,VT>( creator1, creator2 );
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
/*!\brief Macro for the execution of a dense matrix/sparse vector multiplication test case.
*/
#define RUN_DMATSVECMULT_OPERATION_TEST( C1, C2 ) \
   blazetest::mathtest::operations::dmatsvecmult::runTest( C1, C2 )
/*! \endcond */
//*************************************************************************************************

} // namespace dmatsvecmult

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
