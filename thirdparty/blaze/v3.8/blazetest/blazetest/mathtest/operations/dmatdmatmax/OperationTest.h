//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/dmatdmatmax/OperationTest.h
//  \brief Header file for the dense matrix/dense matrix maximum operation test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_DMATDMATMAX_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_DMATDMATMAX_OPERATIONTEST_H_


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
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/Functors.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingScalar.h>
#include <blaze/math/Views.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/constraints/DerivedFrom.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Nor.h>
#include <blaze/util/mpl/Not.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/Random.h>
#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/RemoveCVRef.h>
#include <blazetest/system/LAPACK.h>
#include <blazetest/system/MathTest.h>
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/IsEqual.h>
#include <blazetest/mathtest/MatchAdaptor.h>
#include <blazetest/mathtest/MatchSymmetry.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>


namespace blazetest {

namespace mathtest {

namespace operations {

namespace dmatdmatmax {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the dense matrix/dense matrix maximum operation test.
//
// This class template represents one particular matrix maximum test between two matrices of
// a particular type. The two template arguments \a MT1 and \a MT2 represent the types of the
// left-hand side and right-hand side matrix, respectively.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
class OperationTest
{
 private:
   //**Type definitions****************************************************************************
   using ET1 = blaze::ElementType_t<MT1>;  //!< Element type 1
   using ET2 = blaze::ElementType_t<MT2>;  //!< Element type 2

   using OMT1  = blaze::OppositeType_t<MT1>;    //!< Matrix type 1 with opposite storage order
   using OMT2  = blaze::OppositeType_t<MT2>;    //!< Matrix type 2 with opposite storage order
   using TMT1  = blaze::TransposeType_t<MT1>;   //!< Transpose matrix type 1
   using TMT2  = blaze::TransposeType_t<MT2>;   //!< Transpose matrix type 2
   using TOMT1 = blaze::TransposeType_t<OMT1>;  //!< Transpose matrix type 1 with opposite storage order
   using TOMT2 = blaze::TransposeType_t<OMT2>;  //!< Transpose matrix type 2 with opposite storage order

   //! Dense result type
   using DRE = blaze::MapTrait_t<MT1,MT2,blaze::Max>;

   using DET   = blaze::ElementType_t<DRE>;     //!< Element type of the dense result
   using ODRE  = blaze::OppositeType_t<DRE>;    //!< Dense result type with opposite storage order
   using TDRE  = blaze::TransposeType_t<DRE>;   //!< Transpose dense result type
   using TODRE = blaze::TransposeType_t<ODRE>;  //!< Transpose dense result type with opposite storage order

   //! Sparse result type
   using SRE = MatchAdaptor_t< DRE, blaze::CompressedMatrix<DET,false> >;

   using SET   = blaze::ElementType_t<SRE>;     //!< Element type of the sparse result
   using OSRE  = blaze::OppositeType_t<SRE>;    //!< Sparse result type with opposite storage order
   using TSRE  = blaze::TransposeType_t<SRE>;   //!< Transpose sparse result type
   using TOSRE = blaze::TransposeType_t<OSRE>;  //!< Transpose sparse result type with opposite storage order

   //! Reference result type
   using RT = MatchSymmetry_t< DRE, blaze::DynamicMatrix<blaze::ElementType_t<DRE>,false> >;

   //! Type of the matrix/matrix maximum expression
   using MatMatMaxExprType =
      blaze::RemoveCVRef_t< decltype( max( std::declval<MT1>(), std::declval<MT2>() ) ) >;

   //! Type of the matrix/transpose matrix maximum expression
   using MatTMatMaxExprType =
      blaze::RemoveCVRef_t< decltype( max( std::declval<MT1>(), std::declval<OMT2>() ) ) >;

   //! Type of the transpose matrix/matrix maximum expression
   using TMatMatMaxExprType =
      blaze::RemoveCVRef_t< decltype( max( std::declval<OMT1>(), std::declval<MT2>() ) ) >;

   //! Type of the transpose matrix/transpose matrix maximum expression
   using TMatTMatMaxExprType =
      blaze::RemoveCVRef_t< decltype( max( std::declval<OMT1>(), std::declval<OMT2>() ) ) >;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit OperationTest( const Creator<MT1>& creator1, const Creator<MT2>& creator2 );
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
                          void testInvOperation      ( blaze::TrueType  );
                          void testInvOperation      ( blaze::FalseType );
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
   template< typename LT, typename RT > void checkResults();
   template< typename LT, typename RT > void checkTransposeResults();
   void checkExceptionMessage( const std::exception& ex, const std::string& message );
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
   MT1   lhs_;     //!< The left-hand side dense matrix.
   MT2   rhs_;     //!< The right-hand side dense matrix.
   OMT1  olhs_;    //!< The left-hand side dense matrix with opposite storage order.
   OMT2  orhs_;    //!< The right-hand side dense matrix with opposite storage order.
   DRE   dres_;    //!< The dense result matrix.
   SRE   sres_;    //!< The sparse result matrix.
   ODRE  odres_;   //!< The dense result matrix with opposite storage order.
   OSRE  osres_;   //!< The sparse result matrix with opposite storage order.
   TDRE  tdres_;   //!< The transpose dense result matrix.
   TSRE  tsres_;   //!< The transpose sparse result matrix.
   TODRE todres_;  //!< The transpose dense result matrix with opposite storage order.
   TOSRE tosres_;  //!< The transpose sparse result matrix with opposite storage order.
   RT    ref_;     //!< The reference matrix.
   RT    refres_;  //!< The reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MT1   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MT2   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( OMT1  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( OMT2  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TMT1  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TMT2  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TOMT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TOMT2 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( RT    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( ODRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OSRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TODRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TOSRE );

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( MT1   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( MT2   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OMT1  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OMT2  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TMT1  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TMT2  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TOMT1 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TOMT2 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( RT    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ODRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OSRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TODRE );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TOSRE );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET1, blaze::ElementType_t<OMT1>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET2, blaze::ElementType_t<OMT2>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET1, blaze::ElementType_t<TMT1>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET2, blaze::ElementType_t<TMT2>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET1, blaze::ElementType_t<TOMT1>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET2, blaze::ElementType_t<TOMT2>  );
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
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT1, blaze::OppositeType_t<OMT1>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT2, blaze::OppositeType_t<OMT2>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT1, blaze::TransposeType_t<TMT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT2, blaze::TransposeType_t<TMT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DRE, blaze::OppositeType_t<ODRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DRE, blaze::TransposeType_t<TDRE> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SRE, blaze::OppositeType_t<OSRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SRE, blaze::TransposeType_t<TSRE> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( MatMatMaxExprType, blaze::ResultType_t<MatMatMaxExprType>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MatMatMaxExprType, blaze::OppositeType_t<MatMatMaxExprType>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MatMatMaxExprType, blaze::TransposeType_t<MatMatMaxExprType> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( MatTMatMaxExprType, blaze::ResultType_t<MatTMatMaxExprType>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MatTMatMaxExprType, blaze::OppositeType_t<MatTMatMaxExprType>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MatTMatMaxExprType, blaze::TransposeType_t<MatTMatMaxExprType> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( TMatMatMaxExprType, blaze::ResultType_t<TMatMatMaxExprType>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TMatMatMaxExprType, blaze::OppositeType_t<TMatMatMaxExprType>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TMatMatMaxExprType, blaze::TransposeType_t<TMatMatMaxExprType> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( TMatTMatMaxExprType, blaze::ResultType_t<TMatTMatMaxExprType>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TMatTMatMaxExprType, blaze::OppositeType_t<TMatTMatMaxExprType>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TMatTMatMaxExprType, blaze::TransposeType_t<TMatTMatMaxExprType> );

   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( MatMatMaxExprType  , blaze::BaseType_t<MatMatMaxExprType  > );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( MatTMatMaxExprType , blaze::BaseType_t<MatTMatMaxExprType > );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TMatMatMaxExprType , blaze::BaseType_t<TMatMatMaxExprType > );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TMatTMatMaxExprType, blaze::BaseType_t<TMatTMatMaxExprType> );
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
/*!\brief Constructor for the dense matrix/dense matrix maximum operation test.
//
// \param creator1 The creator for the left-hand side dense matrix of the matrix maximum.
// \param creator2 The creator for the right-hand side dense matrix of the matrix maximum.
// \exception std::runtime_error Operation error detected.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
OperationTest<MT1,MT2>::OperationTest( const Creator<MT1>& creator1, const Creator<MT2>& creator2 )
   : lhs_( creator1() )  // The left-hand side dense matrix
   , rhs_( creator2() )  // The right-hand side dense matrix
   , olhs_( lhs_ )       // The left-hand side dense matrix with opposite storage order
   , orhs_( rhs_ )       // The right-hand side dense matrix with opposite storage order
   , dres_()             // The dense result matrix
   , sres_()             // The sparse result matrix
   , odres_()            // The dense result matrix with opposite storage order
   , osres_()            // The sparse result matrix with opposite storage order
   , tdres_()            // The transpose dense result matrix
   , tsres_()            // The transpose sparse result matrix
   , todres_()           // The transpose dense result matrix with opposite storage order
   , tosres_()           // The transpose sparse result matrix with opposite storage order
   , ref_()              // The reference left-hand side matrix
   , refres_()           // The reference result
   , test_()             // Label of the currently performed test
   , error_()            // Description of the current error type
{
   using namespace blaze;

   using Scalar = blaze::UnderlyingScalar_t<DET>;

   if( lhs_.rows() != rhs_.rows() || lhs_.columns() != rhs_.columns() ) {
      throw std::runtime_error( "Non-matching operands detected" );
   }

   ref_.resize( lhs_.rows(), lhs_.columns() );
   for( size_t i=0UL; i<lhs_.rows(); ++i ) {
      const size_t jbegin( IsUpper<RT>::value ? i : 0UL );
      const size_t jend  ( IsLower<RT>::value ? i+1UL : lhs_.columns() );
      for( size_t j=jbegin; j<jend; ++j ) {
         ref_(i,j) = max( lhs_(i,j), rhs_(i,j) );
      }
   }

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
   testInvOperation( Not_t< IsUniform<DRE> >() );
   testEvalOperation();
   testSerialOperation();
   testNoAliasOperation();
   testNoSIMDOperation();
   testDeclSymOperation( Or_t< IsSquare<DRE>, IsResizable<DRE> >() );
   testDeclHermOperation( Or_t< IsSquare<DRE>, IsResizable<DRE> >() );
   testDeclLowOperation( Or_t< IsSquare<DRE>, IsResizable<DRE> >() );
   testDeclUppOperation( Or_t< IsSquare<DRE>, IsResizable<DRE> >() );
   testDeclDiagOperation( Or_t< IsSquare<DRE>, IsResizable<DRE> >() );
   testSubmatrixOperation( Not_t< IsUniform<DRE> >() );
   testRowOperation( Not_t< IsUniform<DRE> >() );
   testRowsOperation( Nor_t< IsUniform<DRE>, IsSymmetric<DRE>, IsHermitian<DRE> >() );
   testColumnOperation( Not_t< IsUniform<DRE> >() );
   testColumnsOperation( Nor_t< IsUniform<DRE>, IsSymmetric<DRE>, IsHermitian<DRE> >() );
   testBandOperation( Not_t< IsUniform<DRE> >() );
}
//*************************************************************************************************






//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the explicit evaluation.
//
// \return void
// \exception std::runtime_error Evaluation error detected.
//
// This function tests the explicit evaluation. In case any error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testEvaluation()
{
   using blaze::IsRowMajorMatrix;


   //=====================================================================================
   // Testing the evaluation with two row-major matrices
   //=====================================================================================

   {
      const auto res   ( evaluate( max( lhs_, rhs_ ) ) );
      const auto refres( evaluate( ref_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<MT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<MT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
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
      const auto res   ( evaluate( max( eval( lhs_ ), eval( rhs_ ) ) ) );
      const auto refres( evaluate( eval( ref_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<MT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<MT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
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
   // Testing the evaluation with a row-major matrix and a column-major matrix
   //=====================================================================================

   {
      const auto res   ( evaluate( max( lhs_, orhs_ ) ) );
      const auto refres( evaluate( ref_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<MT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<OMT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( orhs_ ).name() << "\n"
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
      const auto res   ( evaluate( max( eval( lhs_ ), eval( orhs_ ) ) ) );
      const auto refres( evaluate( eval( ref_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<MT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<OMT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( orhs_ ).name() << "\n"
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
   // Testing the evaluation with a column-major matrix and a row-major matrix
   //=====================================================================================

   {
      const auto res   ( evaluate( max( olhs_, rhs_ ) ) );
      const auto refres( evaluate( ref_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<OMT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( olhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<MT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
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
      const auto res   ( evaluate( max( eval( olhs_ ), eval( rhs_ ) ) ) );
      const auto refres( evaluate( eval( ref_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<OMT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( olhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<MT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
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
   // Testing the evaluation with two column-major matrices
   //=====================================================================================

   {
      const auto res   ( evaluate( max( olhs_, orhs_ ) ) );
      const auto refres( evaluate( ref_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<OMT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( olhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<OMT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( orhs_ ).name() << "\n"
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
      const auto res   ( evaluate( max( eval( olhs_ ), eval( orhs_ ) ) ) );
      const auto refres( evaluate( eval( ref_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<OMT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( olhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<OMT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
             << "     " << typeid( orhs_ ).name() << "\n"
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
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testElementAccess()
{
   using blaze::equal;


   //=====================================================================================
   // Testing the element access with two row-major matrices
   //=====================================================================================

   if( lhs_.rows() > 0UL && lhs_.columns() > 0UL )
   {
      const size_t m( lhs_.rows()    - 1UL );
      const size_t n( lhs_.columns() - 1UL );

      if( !equal( max( lhs_, rhs_ )(m,n), ref_(m,n) ) ||
          !equal( max( lhs_, rhs_ ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( max( lhs_, eval( rhs_ ) )(m,n), ref_(m,n) ) ||
          !equal( max( lhs_, eval( rhs_ ) ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( max( eval( lhs_ ), rhs_ )(m,n), ref_(m,n) ) ||
          !equal( max( eval( lhs_ ), rhs_ ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( max( eval( lhs_ ), eval( rhs_ ) )(m,n), ref_(m,n) ) ||
          !equal( max( eval( lhs_ ), eval( rhs_ ) ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      max( lhs_, rhs_ ).at( 0UL, lhs_.columns() );

      std::ostringstream oss;
      oss << " Test : Checked element access of maximum expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side row-major dense matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Right-hand side row-major dense matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      max( lhs_, rhs_ ).at( lhs_.rows(), 0UL );

      std::ostringstream oss;
      oss << " Test : Checked element access of maximum expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side row-major dense matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Right-hand side row-major dense matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}


   //=====================================================================================
   // Testing the element access with a row-major matrix and a column-major matrix
   //=====================================================================================

   if( lhs_.rows() > 0UL && lhs_.columns() > 0UL )
   {
      const size_t m( lhs_.rows()    - 1UL );
      const size_t n( lhs_.columns() - 1UL );

      if( !equal( max( lhs_, orhs_ )(m,n), ref_(m,n) ) ||
          !equal( max( lhs_, orhs_ ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( max( lhs_, eval( orhs_ ) )(m,n), ref_(m,n) ) ||
          !equal( max( lhs_, eval( orhs_ ) ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( max( eval( lhs_ ), orhs_ )(m,n), ref_(m,n) ) ||
          !equal( max( eval( lhs_ ), orhs_ ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( max( eval( lhs_ ), eval( orhs_ ) )(m,n), ref_(m,n) ) ||
          !equal( max( eval( lhs_ ), eval( orhs_ ) ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      max( lhs_, orhs_ ).at( 0UL, lhs_.columns() );

      std::ostringstream oss;
      oss << " Test : Checked element access of maximum expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side row-major dense matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Right-hand side column-major dense matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      max( lhs_, orhs_ ).at( lhs_.rows(), 0UL );

      std::ostringstream oss;
      oss << " Test : Checked element access of maximum expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side row-major dense matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Right-hand side column-major dense matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}


   //=====================================================================================
   // Testing the element access with a column-major matrix and a row-major matrix
   //=====================================================================================

   if( olhs_.rows() > 0UL && olhs_.columns() > 0UL )
   {
      const size_t m( olhs_.rows()    - 1UL );
      const size_t n( olhs_.columns() - 1UL );

      if( !equal( max( olhs_, rhs_ )(m,n), ref_(m,n) ) ||
          !equal( max( olhs_, rhs_ ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( max( olhs_, eval( rhs_ ) )(m,n), ref_(m,n) ) ||
          !equal( max( olhs_, eval( rhs_ ) ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( max( eval( olhs_ ), rhs_ )(m,n), ref_(m,n) ) ||
          !equal( max( eval( olhs_ ), rhs_ ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( max( eval( olhs_ ), eval( rhs_ ) )(m,n), ref_(m,n) ) ||
          !equal( max( eval( olhs_ ), eval( rhs_ ) ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side row-major dense matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      max( olhs_, rhs_ ).at( 0UL, lhs_.columns() );

      std::ostringstream oss;
      oss << " Test : Checked element access of maximum expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side column-major dense matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Right-hand side row-major dense matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      max( olhs_, rhs_ ).at( lhs_.rows(), 0UL );

      std::ostringstream oss;
      oss << " Test : Checked element access of maximum expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side column-major dense matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Right-hand side row-major dense matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}


   //=====================================================================================
   // Testing the element access with two column-major matrices
   //=====================================================================================

   if( olhs_.rows() > 0UL && olhs_.columns() > 0UL )
   {
      const size_t m( olhs_.rows()    - 1UL );
      const size_t n( olhs_.columns() - 1UL );

      if( !equal( max( olhs_, orhs_ )(m,n), ref_(m,n) ) ||
          !equal( max( olhs_, orhs_ ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( max( olhs_, eval( orhs_ ) )(m,n), ref_(m,n) ) ||
          !equal( max( olhs_, eval( orhs_ ) ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( max( eval( olhs_ ), orhs_ )(m,n), ref_(m,n) ) ||
          !equal( max( eval( olhs_ ), orhs_ ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( max( eval( olhs_ ), eval( orhs_ ) )(m,n), ref_(m,n) ) ||
          !equal( max( eval( olhs_ ), eval( orhs_ ) ).at(m,n), ref_.at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated maximum expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side column-major dense matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      max( olhs_, orhs_ ).at( 0UL, lhs_.columns() );

      std::ostringstream oss;
      oss << " Test : Checked element access of maximum expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side column-major dense matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Right-hand side column-major dense matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      max( olhs_, orhs_ ).at( lhs_.rows(), 0UL );

      std::ostringstream oss;
      oss << " Test : Checked element access of maximum expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side column-major dense matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Right-hand side column-major dense matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from
// the maximum operation or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Maximum
      //=====================================================================================

      // Maximum with the given matrices
      {
         test_  = "Maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = max( lhs_, rhs_ );
            odres_  = max( lhs_, rhs_ );
            sres_   = max( lhs_, rhs_ );
            osres_  = max( lhs_, rhs_ );
            refres_ = ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = max( lhs_, orhs_ );
            odres_  = max( lhs_, orhs_ );
            sres_   = max( lhs_, orhs_ );
            osres_  = max( lhs_, orhs_ );
            refres_ = ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = max( olhs_, rhs_ );
            odres_  = max( olhs_, rhs_ );
            sres_   = max( olhs_, rhs_ );
            osres_  = max( olhs_, rhs_ );
            refres_ = ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = max( olhs_, orhs_ );
            odres_  = max( olhs_, orhs_ );
            sres_   = max( olhs_, orhs_ );
            osres_  = max( olhs_, orhs_ );
            refres_ = ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Maximum with evaluated matrices
      {
         test_  = "Maximum with evaluated matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = max( eval( lhs_ ), eval( rhs_ ) );
            odres_  = max( eval( lhs_ ), eval( rhs_ ) );
            sres_   = max( eval( lhs_ ), eval( rhs_ ) );
            osres_  = max( eval( lhs_ ), eval( rhs_ ) );
            refres_ = eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = max( eval( lhs_ ), eval( orhs_ ) );
            odres_  = max( eval( lhs_ ), eval( orhs_ ) );
            sres_   = max( eval( lhs_ ), eval( orhs_ ) );
            osres_  = max( eval( lhs_ ), eval( orhs_ ) );
            refres_ = eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = max( eval( olhs_ ), eval( rhs_ ) );
            odres_  = max( eval( olhs_ ), eval( rhs_ ) );
            sres_   = max( eval( olhs_ ), eval( rhs_ ) );
            osres_  = max( eval( olhs_ ), eval( rhs_ ) );
            refres_ = eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = max( eval( olhs_ ), eval( orhs_ ) );
            odres_  = max( eval( olhs_ ), eval( orhs_ ) );
            sres_   = max( eval( olhs_ ), eval( orhs_ ) );
            osres_  = max( eval( olhs_ ), eval( orhs_ ) );
            refres_ = eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Maximum with addition assignment
      //=====================================================================================

      // Maximum with addition assignment with the given matrices
      {
         test_  = "Maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += max( lhs_, rhs_ );
            odres_  += max( lhs_, rhs_ );
            sres_   += max( lhs_, rhs_ );
            osres_  += max( lhs_, rhs_ );
            refres_ += ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += max( lhs_, orhs_ );
            odres_  += max( lhs_, orhs_ );
            sres_   += max( lhs_, orhs_ );
            osres_  += max( lhs_, orhs_ );
            refres_ += ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += max( olhs_, rhs_ );
            odres_  += max( olhs_, rhs_ );
            sres_   += max( olhs_, rhs_ );
            osres_  += max( olhs_, rhs_ );
            refres_ += ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += max( olhs_, orhs_ );
            odres_  += max( olhs_, orhs_ );
            sres_   += max( olhs_, orhs_ );
            osres_  += max( olhs_, orhs_ );
            refres_ += ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Maximum with addition assignment with evaluated matrices
      {
         test_  = "Maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += max( eval( lhs_ ), eval( rhs_ ) );
            odres_  += max( eval( lhs_ ), eval( rhs_ ) );
            sres_   += max( eval( lhs_ ), eval( rhs_ ) );
            osres_  += max( eval( lhs_ ), eval( rhs_ ) );
            refres_ += eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += max( eval( lhs_ ), eval( orhs_ ) );
            odres_  += max( eval( lhs_ ), eval( orhs_ ) );
            sres_   += max( eval( lhs_ ), eval( orhs_ ) );
            osres_  += max( eval( lhs_ ), eval( orhs_ ) );
            refres_ += eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += max( eval( olhs_ ), eval( rhs_ ) );
            odres_  += max( eval( olhs_ ), eval( rhs_ ) );
            sres_   += max( eval( olhs_ ), eval( rhs_ ) );
            osres_  += max( eval( olhs_ ), eval( rhs_ ) );
            refres_ += eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += max( eval( olhs_ ), eval( orhs_ ) );
            odres_  += max( eval( olhs_ ), eval( orhs_ ) );
            sres_   += max( eval( olhs_ ), eval( orhs_ ) );
            osres_  += max( eval( olhs_ ), eval( orhs_ ) );
            refres_ += eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Maximum with subtraction assignment with the given matrices
      //=====================================================================================

      // Maximum with subtraction assignment with the given matrices
      {
         test_  = "Maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= max( lhs_, rhs_ );
            odres_  -= max( lhs_, rhs_ );
            sres_   -= max( lhs_, rhs_ );
            osres_  -= max( lhs_, rhs_ );
            refres_ -= ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= max( lhs_, orhs_ );
            odres_  -= max( lhs_, orhs_ );
            sres_   -= max( lhs_, orhs_ );
            osres_  -= max( lhs_, orhs_ );
            refres_ -= ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= max( olhs_, rhs_ );
            odres_  -= max( olhs_, rhs_ );
            sres_   -= max( olhs_, rhs_ );
            osres_  -= max( olhs_, rhs_ );
            refres_ -= ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= max( olhs_, orhs_ );
            odres_  -= max( olhs_, orhs_ );
            sres_   -= max( olhs_, orhs_ );
            osres_  -= max( olhs_, orhs_ );
            refres_ -= ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= max( eval( lhs_ ), eval( rhs_ ) );
            odres_  -= max( eval( lhs_ ), eval( rhs_ ) );
            sres_   -= max( eval( lhs_ ), eval( rhs_ ) );
            osres_  -= max( eval( lhs_ ), eval( rhs_ ) );
            refres_ -= eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= max( eval( lhs_ ), eval( orhs_ ) );
            odres_  -= max( eval( lhs_ ), eval( orhs_ ) );
            sres_   -= max( eval( lhs_ ), eval( orhs_ ) );
            osres_  -= max( eval( lhs_ ), eval( orhs_ ) );
            refres_ -= eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= max( eval( olhs_ ), eval( rhs_ ) );
            odres_  -= max( eval( olhs_ ), eval( rhs_ ) );
            sres_   -= max( eval( olhs_ ), eval( rhs_ ) );
            osres_  -= max( eval( olhs_ ), eval( rhs_ ) );
            refres_ -= eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= max( eval( olhs_ ), eval( orhs_ ) );
            odres_  -= max( eval( olhs_ ), eval( orhs_ ) );
            sres_   -= max( eval( olhs_ ), eval( orhs_ ) );
            osres_  -= max( eval( olhs_ ), eval( orhs_ ) );
            refres_ -= eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Maximum with Schur product assignment
      //=====================================================================================

      // Maximum with Schur product assignment with the given matrices
      {
         test_  = "Maximum with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= max( lhs_, rhs_ );
            odres_  %= max( lhs_, rhs_ );
            sres_   %= max( lhs_, rhs_ );
            osres_  %= max( lhs_, rhs_ );
            refres_ %= ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= max( lhs_, orhs_ );
            odres_  %= max( lhs_, orhs_ );
            sres_   %= max( lhs_, orhs_ );
            osres_  %= max( lhs_, orhs_ );
            refres_ %= ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= max( olhs_, rhs_ );
            odres_  %= max( olhs_, rhs_ );
            sres_   %= max( olhs_, rhs_ );
            osres_  %= max( olhs_, rhs_ );
            refres_ %= ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= max( olhs_, orhs_ );
            odres_  %= max( olhs_, orhs_ );
            sres_   %= max( olhs_, orhs_ );
            osres_  %= max( olhs_, orhs_ );
            refres_ %= ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Maximum with Schur product assignment with evaluated matrices
      {
         test_  = "Maximum with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= max( eval( lhs_ ), eval( rhs_ ) );
            odres_  %= max( eval( lhs_ ), eval( rhs_ ) );
            sres_   %= max( eval( lhs_ ), eval( rhs_ ) );
            osres_  %= max( eval( lhs_ ), eval( rhs_ ) );
            refres_ %= eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= max( eval( lhs_ ), eval( orhs_ ) );
            odres_  %= max( eval( lhs_ ), eval( orhs_ ) );
            sres_   %= max( eval( lhs_ ), eval( orhs_ ) );
            osres_  %= max( eval( lhs_ ), eval( orhs_ ) );
            refres_ %= eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= max( eval( olhs_ ), eval( rhs_ ) );
            odres_  %= max( eval( olhs_ ), eval( rhs_ ) );
            sres_   %= max( eval( olhs_ ), eval( rhs_ ) );
            osres_  %= max( eval( olhs_ ), eval( rhs_ ) );
            refres_ %= eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= max( eval( olhs_ ), eval( orhs_ ) );
            odres_  %= max( eval( olhs_ ), eval( orhs_ ) );
            sres_   %= max( eval( olhs_ ), eval( orhs_ ) );
            osres_  %= max( eval( olhs_ ), eval( orhs_ ) );
            refres_ %= eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the negated dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the negated matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      //=====================================================================================
      // Negated maximum
      //=====================================================================================

      // Negated maximum with the given matrices
      {
         test_  = "Negated maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = -max( lhs_, rhs_ );
            odres_  = -max( lhs_, rhs_ );
            sres_   = -max( lhs_, rhs_ );
            osres_  = -max( lhs_, rhs_ );
            refres_ = -ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = -max( lhs_, orhs_ );
            odres_  = -max( lhs_, orhs_ );
            sres_   = -max( lhs_, orhs_ );
            osres_  = -max( lhs_, orhs_ );
            refres_ = -ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = -max( olhs_, rhs_ );
            odres_  = -max( olhs_, rhs_ );
            sres_   = -max( olhs_, rhs_ );
            osres_  = -max( olhs_, rhs_ );
            refres_ = -ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = -max( olhs_, orhs_ );
            odres_  = -max( olhs_, orhs_ );
            sres_   = -max( olhs_, orhs_ );
            osres_  = -max( olhs_, orhs_ );
            refres_ = -ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Negated maximum with evaluated matrices
      {
         test_  = "Negated maximum with evaluated matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = -max( eval( lhs_ ), eval( rhs_ ) );
            odres_  = -max( eval( lhs_ ), eval( rhs_ ) );
            sres_   = -max( eval( lhs_ ), eval( rhs_ ) );
            osres_  = -max( eval( lhs_ ), eval( rhs_ ) );
            refres_ = -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = -max( eval( lhs_ ), eval( orhs_ ) );
            odres_  = -max( eval( lhs_ ), eval( orhs_ ) );
            sres_   = -max( eval( lhs_ ), eval( orhs_ ) );
            osres_  = -max( eval( lhs_ ), eval( orhs_ ) );
            refres_ = -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = -max( eval( olhs_ ), eval( rhs_ ) );
            odres_  = -max( eval( olhs_ ), eval( rhs_ ) );
            sres_   = -max( eval( olhs_ ), eval( rhs_ ) );
            osres_  = -max( eval( olhs_ ), eval( rhs_ ) );
            refres_ = -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = -max( eval( olhs_ ), eval( orhs_ ) );
            odres_  = -max( eval( olhs_ ), eval( orhs_ ) );
            sres_   = -max( eval( olhs_ ), eval( orhs_ ) );
            osres_  = -max( eval( olhs_ ), eval( orhs_ ) );
            refres_ = -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Negated maximum with addition assignment
      //=====================================================================================

      // Negated maximum with addition assignment with the given matrices
      {
         test_  = "Negated maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -max( lhs_, rhs_ );
            odres_  += -max( lhs_, rhs_ );
            sres_   += -max( lhs_, rhs_ );
            osres_  += -max( lhs_, rhs_ );
            refres_ += -ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += -max( lhs_, orhs_ );
            odres_  += -max( lhs_, orhs_ );
            sres_   += -max( lhs_, orhs_ );
            osres_  += -max( lhs_, orhs_ );
            refres_ += -ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += -max( olhs_, rhs_ );
            odres_  += -max( olhs_, rhs_ );
            sres_   += -max( olhs_, rhs_ );
            osres_  += -max( olhs_, rhs_ );
            refres_ += -ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += -max( olhs_, orhs_ );
            odres_  += -max( olhs_, orhs_ );
            sres_   += -max( olhs_, orhs_ );
            osres_  += -max( olhs_, orhs_ );
            refres_ += -ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Negated maximum with addition assignment with the given matrices
      {
         test_  = "Negated maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -max( eval( lhs_ ), eval( rhs_ ) );
            odres_  += -max( eval( lhs_ ), eval( rhs_ ) );
            sres_   += -max( eval( lhs_ ), eval( rhs_ ) );
            osres_  += -max( eval( lhs_ ), eval( rhs_ ) );
            refres_ += -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += -max( eval( lhs_ ), eval( orhs_ ) );
            odres_  += -max( eval( lhs_ ), eval( orhs_ ) );
            sres_   += -max( eval( lhs_ ), eval( orhs_ ) );
            osres_  += -max( eval( lhs_ ), eval( orhs_ ) );
            refres_ += -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += -max( eval( olhs_ ), eval( rhs_ ) );
            odres_  += -max( eval( olhs_ ), eval( rhs_ ) );
            sres_   += -max( eval( olhs_ ), eval( rhs_ ) );
            osres_  += -max( eval( olhs_ ), eval( rhs_ ) );
            refres_ += -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += -max( eval( olhs_ ), eval( orhs_ ) );
            odres_  += -max( eval( olhs_ ), eval( orhs_ ) );
            sres_   += -max( eval( olhs_ ), eval( orhs_ ) );
            osres_  += -max( eval( olhs_ ), eval( orhs_ ) );
            refres_ += -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Negated maximum with subtraction assignment
      //=====================================================================================

      // Negated maximum with subtraction assignment with the given matrices
      {
         test_  = "Negated maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -max( lhs_, rhs_ );
            odres_  -= -max( lhs_, rhs_ );
            sres_   -= -max( lhs_, rhs_ );
            osres_  -= -max( lhs_, rhs_ );
            refres_ -= -ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= -max( lhs_, orhs_ );
            odres_  -= -max( lhs_, orhs_ );
            sres_   -= -max( lhs_, orhs_ );
            osres_  -= -max( lhs_, orhs_ );
            refres_ -= -ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= -max( olhs_, rhs_ );
            odres_  -= -max( olhs_, rhs_ );
            sres_   -= -max( olhs_, rhs_ );
            osres_  -= -max( olhs_, rhs_ );
            refres_ -= -ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= -max( olhs_, orhs_ );
            odres_  -= -max( olhs_, orhs_ );
            sres_   -= -max( olhs_, orhs_ );
            osres_  -= -max( olhs_, orhs_ );
            refres_ -= -ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Negated maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Negated maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -max( eval( lhs_ ), eval( rhs_ ) );
            odres_  -= -max( eval( lhs_ ), eval( rhs_ ) );
            sres_   -= -max( eval( lhs_ ), eval( rhs_ ) );
            osres_  -= -max( eval( lhs_ ), eval( rhs_ ) );
            refres_ -= -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= -max( eval( lhs_ ), eval( orhs_ ) );
            odres_  -= -max( eval( lhs_ ), eval( orhs_ ) );
            sres_   -= -max( eval( lhs_ ), eval( orhs_ ) );
            osres_  -= -max( eval( lhs_ ), eval( orhs_ ) );
            refres_ -= -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= -max( eval( olhs_ ), eval( rhs_ ) );
            odres_  -= -max( eval( olhs_ ), eval( rhs_ ) );
            sres_   -= -max( eval( olhs_ ), eval( rhs_ ) );
            osres_  -= -max( eval( olhs_ ), eval( rhs_ ) );
            refres_ -= -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= -max( eval( olhs_ ), eval( orhs_ ) );
            odres_  -= -max( eval( olhs_ ), eval( orhs_ ) );
            sres_   -= -max( eval( olhs_ ), eval( orhs_ ) );
            osres_  -= -max( eval( olhs_ ), eval( orhs_ ) );
            refres_ -= -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Negated maximum with Schur product assignment
      //=====================================================================================

      // Negated maximum with Schur product assignment with the given matrices
      {
         test_  = "Negated maximum with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= -max( lhs_, rhs_ );
            odres_  %= -max( lhs_, rhs_ );
            sres_   %= -max( lhs_, rhs_ );
            osres_  %= -max( lhs_, rhs_ );
            refres_ %= -ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= -max( lhs_, orhs_ );
            odres_  %= -max( lhs_, orhs_ );
            sres_   %= -max( lhs_, orhs_ );
            osres_  %= -max( lhs_, orhs_ );
            refres_ %= -ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= -max( olhs_, rhs_ );
            odres_  %= -max( olhs_, rhs_ );
            sres_   %= -max( olhs_, rhs_ );
            osres_  %= -max( olhs_, rhs_ );
            refres_ %= -ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= -max( olhs_, orhs_ );
            odres_  %= -max( olhs_, orhs_ );
            sres_   %= -max( olhs_, orhs_ );
            osres_  %= -max( olhs_, orhs_ );
            refres_ %= -ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Negated maximum with Schur product assignment with the given matrices
      {
         test_  = "Negated maximum with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= -max( eval( lhs_ ), eval( rhs_ ) );
            odres_  %= -max( eval( lhs_ ), eval( rhs_ ) );
            sres_   %= -max( eval( lhs_ ), eval( rhs_ ) );
            osres_  %= -max( eval( lhs_ ), eval( rhs_ ) );
            refres_ %= -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= -max( eval( lhs_ ), eval( orhs_ ) );
            odres_  %= -max( eval( lhs_ ), eval( orhs_ ) );
            sres_   %= -max( eval( lhs_ ), eval( orhs_ ) );
            osres_  %= -max( eval( lhs_ ), eval( orhs_ ) );
            refres_ %= -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= -max( eval( olhs_ ), eval( rhs_ ) );
            odres_  %= -max( eval( olhs_ ), eval( rhs_ ) );
            sres_   %= -max( eval( olhs_ ), eval( rhs_ ) );
            osres_  %= -max( eval( olhs_ ), eval( rhs_ ) );
            refres_ %= -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= -max( eval( olhs_ ), eval( orhs_ ) );
            odres_  %= -max( eval( olhs_ ), eval( orhs_ ) );
            sres_   %= -max( eval( olhs_ ), eval( orhs_ ) );
            osres_  %= -max( eval( olhs_ ), eval( orhs_ ) );
            refres_ %= -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the scaled dense matrix/dense matrix maximum operation.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the scaled matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
template< typename T >    // Type of the scalar
void OperationTest<MT1,MT2>::testScaledOperation( T scalar )
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
            dres_   = max( lhs_, rhs_ );
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

         checkResults<MT1,MT2>();
      }


      //=====================================================================================
      // Self-scaling (M=M*s)
      //=====================================================================================

      // Self-scaling (M=M*s)
      {
         test_ = "Self-scaling (M=M*s)";

         try {
            dres_   = max( lhs_, rhs_ );
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

         checkResults<MT1,MT2>();
      }


      //=====================================================================================
      // Self-scaling (M=s*M)
      //=====================================================================================

      // Self-scaling (M=s*M)
      {
         test_ = "Self-scaling (M=s*M)";

         try {
            dres_   = max( lhs_, rhs_ );
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

         checkResults<MT1,MT2>();
      }


      //=====================================================================================
      // Self-scaling (M/=s)
      //=====================================================================================

      // Self-scaling (M/=s)
      {
         test_ = "Self-scaling (M/=s)";

         try {
            dres_   = max( lhs_, rhs_ );
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

         checkResults<MT1,MT2>();
      }


      //=====================================================================================
      // Self-scaling (M=M/s)
      //=====================================================================================

      // Self-scaling (M=M/s)
      {
         test_ = "Self-scaling (M=M/s)";

         try {
            dres_   = max( lhs_, rhs_ );
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

         checkResults<MT1,MT2>();
      }


      //=====================================================================================
      // Scaled maximum (s*OP)
      //=====================================================================================

      // Scaled maximum with the given matrices
      {
         test_  = "Scaled maximum with the given matrices (s*OP)";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = scalar * max( lhs_, rhs_ );
            odres_  = scalar * max( lhs_, rhs_ );
            sres_   = scalar * max( lhs_, rhs_ );
            osres_  = scalar * max( lhs_, rhs_ );
            refres_ = scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = scalar * max( lhs_, orhs_ );
            odres_  = scalar * max( lhs_, orhs_ );
            sres_   = scalar * max( lhs_, orhs_ );
            osres_  = scalar * max( lhs_, orhs_ );
            refres_ = scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = scalar * max( olhs_, rhs_ );
            odres_  = scalar * max( olhs_, rhs_ );
            sres_   = scalar * max( olhs_, rhs_ );
            osres_  = scalar * max( olhs_, rhs_ );
            refres_ = scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = scalar * max( olhs_, orhs_ );
            odres_  = scalar * max( olhs_, orhs_ );
            sres_   = scalar * max( olhs_, orhs_ );
            osres_  = scalar * max( olhs_, orhs_ );
            refres_ = scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled maximum with evaluated matrices
      {
         test_  = "Scaled maximum with evaluated matrices (s*OP)";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = scalar * max( eval( lhs_ ), eval( rhs_ ) );
            odres_  = scalar * max( eval( lhs_ ), eval( rhs_ ) );
            sres_   = scalar * max( eval( lhs_ ), eval( rhs_ ) );
            osres_  = scalar * max( eval( lhs_ ), eval( rhs_ ) );
            refres_ = scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = scalar * max( eval( lhs_ ), eval( orhs_ ) );
            odres_  = scalar * max( eval( lhs_ ), eval( orhs_ ) );
            sres_   = scalar * max( eval( lhs_ ), eval( orhs_ ) );
            osres_  = scalar * max( eval( lhs_ ), eval( orhs_ ) );
            refres_ = scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = scalar * max( eval( olhs_ ), eval( rhs_ ) );
            odres_  = scalar * max( eval( olhs_ ), eval( rhs_ ) );
            sres_   = scalar * max( eval( olhs_ ), eval( rhs_ ) );
            osres_  = scalar * max( eval( olhs_ ), eval( rhs_ ) );
            refres_ = scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = scalar * max( eval( olhs_ ), eval( orhs_ ) );
            odres_  = scalar * max( eval( olhs_ ), eval( orhs_ ) );
            sres_   = scalar * max( eval( olhs_ ), eval( orhs_ ) );
            osres_  = scalar * max( eval( olhs_ ), eval( orhs_ ) );
            refres_ = scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled maximum (OP*s)
      //=====================================================================================

      // Scaled maximum with the given matrices
      {
         test_  = "Scaled maximum with the given matrices (OP*s)";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = max( lhs_, rhs_ ) * scalar;
            odres_  = max( lhs_, rhs_ ) * scalar;
            sres_   = max( lhs_, rhs_ ) * scalar;
            osres_  = max( lhs_, rhs_ ) * scalar;
            refres_ = ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = max( lhs_, orhs_ ) * scalar;
            odres_  = max( lhs_, orhs_ ) * scalar;
            sres_   = max( lhs_, orhs_ ) * scalar;
            osres_  = max( lhs_, orhs_ ) * scalar;
            refres_ = ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = max( olhs_, rhs_ ) * scalar;
            odres_  = max( olhs_, rhs_ ) * scalar;
            sres_   = max( olhs_, rhs_ ) * scalar;
            osres_  = max( olhs_, rhs_ ) * scalar;
            refres_ = ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = max( olhs_, orhs_ ) * scalar;
            odres_  = max( olhs_, orhs_ ) * scalar;
            sres_   = max( olhs_, orhs_ ) * scalar;
            osres_  = max( olhs_, orhs_ ) * scalar;
            refres_ = ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled maximum with evaluated matrices
      {
         test_  = "Scaled maximum with evaluated matrices (OP*s)";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            odres_  = max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            sres_   = max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            osres_  = max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            refres_ = eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            odres_  = max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            sres_   = max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            osres_  = max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            refres_ = eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            odres_  = max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            sres_   = max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            osres_  = max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            refres_ = eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            odres_  = max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            sres_   = max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            osres_  = max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            refres_ = eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled maximum (OP/s)
      //=====================================================================================

      // Scaled maximum with the given matrices
      {
         test_  = "Scaled maximum with the given matrices (OP/s)";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = max( lhs_, rhs_ ) / scalar;
            odres_  = max( lhs_, rhs_ ) / scalar;
            sres_   = max( lhs_, rhs_ ) / scalar;
            osres_  = max( lhs_, rhs_ ) / scalar;
            refres_ = ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = max( lhs_, orhs_ ) / scalar;
            odres_  = max( lhs_, orhs_ ) / scalar;
            sres_   = max( lhs_, orhs_ ) / scalar;
            osres_  = max( lhs_, orhs_ ) / scalar;
            refres_ = ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = max( olhs_, rhs_ ) / scalar;
            odres_  = max( olhs_, rhs_ ) / scalar;
            sres_   = max( olhs_, rhs_ ) / scalar;
            osres_  = max( olhs_, rhs_ ) / scalar;
            refres_ = ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = max( olhs_, orhs_ ) / scalar;
            odres_  = max( olhs_, orhs_ ) / scalar;
            sres_   = max( olhs_, orhs_ ) / scalar;
            osres_  = max( olhs_, orhs_ ) / scalar;
            refres_ = ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled maximum with evaluated matrices
      {
         test_  = "Scaled maximum with evaluated matrices (OP/s)";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            odres_  = max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            sres_   = max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            osres_  = max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            refres_ = eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            odres_  = max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            sres_   = max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            osres_  = max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            refres_ = eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            odres_  = max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            sres_   = max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            osres_  = max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            refres_ = eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            odres_  = max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            sres_   = max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            osres_  = max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            refres_ = eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled maximum with addition assignment (s*OP)
      //=====================================================================================

      // Scaled maximum with addition assignment with the given matrices
      {
         test_  = "Scaled maximum with addition assignment with the given matrices (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * max( lhs_, rhs_ );
            odres_  += scalar * max( lhs_, rhs_ );
            sres_   += scalar * max( lhs_, rhs_ );
            osres_  += scalar * max( lhs_, rhs_ );
            refres_ += scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += scalar * max( lhs_, orhs_ );
            odres_  += scalar * max( lhs_, orhs_ );
            sres_   += scalar * max( lhs_, orhs_ );
            osres_  += scalar * max( lhs_, orhs_ );
            refres_ += scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += scalar * max( olhs_, rhs_ );
            odres_  += scalar * max( olhs_, rhs_ );
            sres_   += scalar * max( olhs_, rhs_ );
            osres_  += scalar * max( olhs_, rhs_ );
            refres_ += scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += scalar * max( olhs_, orhs_ );
            odres_  += scalar * max( olhs_, orhs_ );
            sres_   += scalar * max( olhs_, orhs_ );
            osres_  += scalar * max( olhs_, orhs_ );
            refres_ += scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled maximum with addition assignment with evaluated matrices
      {
         test_  = "Scaled maximum with addition assignment with evaluated matrices (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * max( eval( lhs_ ), eval( rhs_ ) );
            odres_  += scalar * max( eval( lhs_ ), eval( rhs_ ) );
            sres_   += scalar * max( eval( lhs_ ), eval( rhs_ ) );
            osres_  += scalar * max( eval( lhs_ ), eval( rhs_ ) );
            refres_ += scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += scalar * max( eval( lhs_ ), eval( orhs_ ) );
            odres_  += scalar * max( eval( lhs_ ), eval( orhs_ ) );
            sres_   += scalar * max( eval( lhs_ ), eval( orhs_ ) );
            osres_  += scalar * max( eval( lhs_ ), eval( orhs_ ) );
            refres_ += scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += scalar * max( eval( olhs_ ), eval( rhs_ ) );
            odres_  += scalar * max( eval( olhs_ ), eval( rhs_ ) );
            sres_   += scalar * max( eval( olhs_ ), eval( rhs_ ) );
            osres_  += scalar * max( eval( olhs_ ), eval( rhs_ ) );
            refres_ += scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += scalar * max( eval( olhs_ ), eval( orhs_ ) );
            odres_  += scalar * max( eval( olhs_ ), eval( orhs_ ) );
            sres_   += scalar * max( eval( olhs_ ), eval( orhs_ ) );
            osres_  += scalar * max( eval( olhs_ ), eval( orhs_ ) );
            refres_ += scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled maximum with addition assignment (OP*s)
      //=====================================================================================

      // Scaled maximum with addition assignment with the given matrices
      {
         test_  = "Scaled maximum with addition assignment with the given matrices (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += max( lhs_, rhs_ ) * scalar;
            odres_  += max( lhs_, rhs_ ) * scalar;
            sres_   += max( lhs_, rhs_ ) * scalar;
            osres_  += max( lhs_, rhs_ ) * scalar;
            refres_ += ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += max( lhs_, orhs_ ) * scalar;
            odres_  += max( lhs_, orhs_ ) * scalar;
            sres_   += max( lhs_, orhs_ ) * scalar;
            osres_  += max( lhs_, orhs_ ) * scalar;
            refres_ += ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += max( olhs_, rhs_ ) * scalar;
            odres_  += max( olhs_, rhs_ ) * scalar;
            sres_   += max( olhs_, rhs_ ) * scalar;
            osres_  += max( olhs_, rhs_ ) * scalar;
            refres_ += ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += max( olhs_, orhs_ ) * scalar;
            odres_  += max( olhs_, orhs_ ) * scalar;
            sres_   += max( olhs_, orhs_ ) * scalar;
            osres_  += max( olhs_, orhs_ ) * scalar;
            refres_ += ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled maximum with addition assignment with evaluated matrices
      {
         test_  = "Scaled maximum with addition assignment with evaluated matrices (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            odres_  += max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            sres_   += max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            osres_  += max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            refres_ += eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            odres_  += max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            sres_   += max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            osres_  += max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            refres_ += eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            odres_  += max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            sres_   += max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            osres_  += max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            refres_ += eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            odres_  += max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            sres_   += max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            osres_  += max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            refres_ += eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled maximum with addition assignment (OP/s)
      //=====================================================================================

      // Scaled maximum with addition assignment with the given matrices
      {
         test_  = "Scaled maximum with addition assignment with the given matrices (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += max( lhs_, rhs_ ) / scalar;
            odres_  += max( lhs_, rhs_ ) / scalar;
            sres_   += max( lhs_, rhs_ ) / scalar;
            osres_  += max( lhs_, rhs_ ) / scalar;
            refres_ += ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += max( lhs_, orhs_ ) / scalar;
            odres_  += max( lhs_, orhs_ ) / scalar;
            sres_   += max( lhs_, orhs_ ) / scalar;
            osres_  += max( lhs_, orhs_ ) / scalar;
            refres_ += ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += max( olhs_, rhs_ ) / scalar;
            odres_  += max( olhs_, rhs_ ) / scalar;
            sres_   += max( olhs_, rhs_ ) / scalar;
            osres_  += max( olhs_, rhs_ ) / scalar;
            refres_ += ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += max( olhs_, orhs_ ) / scalar;
            odres_  += max( olhs_, orhs_ ) / scalar;
            sres_   += max( olhs_, orhs_ ) / scalar;
            osres_  += max( olhs_, orhs_ ) / scalar;
            refres_ += ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled maximum with addition assignment with evaluated matrices
      {
         test_  = "Scaled maximum with addition assignment with evaluated matrices (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            odres_  += max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            sres_   += max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            osres_  += max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            refres_ += eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            odres_  += max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            sres_   += max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            osres_  += max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            refres_ += eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            odres_  += max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            sres_   += max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            osres_  += max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            refres_ += eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            odres_  += max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            sres_   += max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            osres_  += max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            refres_ += eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled maximum with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled maximum with subtraction assignment with the given matrices
      {
         test_  = "Scaled maximum with subtraction assignment with the given matrices (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * max( lhs_, rhs_ );
            odres_  -= scalar * max( lhs_, rhs_ );
            sres_   -= scalar * max( lhs_, rhs_ );
            osres_  -= scalar * max( lhs_, rhs_ );
            refres_ -= scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= scalar * max( lhs_, orhs_ );
            odres_  -= scalar * max( lhs_, orhs_ );
            sres_   -= scalar * max( lhs_, orhs_ );
            osres_  -= scalar * max( lhs_, orhs_ );
            refres_ -= scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= scalar * max( olhs_, rhs_ );
            odres_  -= scalar * max( olhs_, rhs_ );
            sres_   -= scalar * max( olhs_, rhs_ );
            osres_  -= scalar * max( olhs_, rhs_ );
            refres_ -= scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= scalar * max( olhs_, orhs_ );
            odres_  -= scalar * max( olhs_, orhs_ );
            sres_   -= scalar * max( olhs_, orhs_ );
            osres_  -= scalar * max( olhs_, orhs_ );
            refres_ -= scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Scaled maximum with subtraction assignment with evaluated matrices (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * max( eval( lhs_ ), eval( rhs_ ) );
            odres_  -= scalar * max( eval( lhs_ ), eval( rhs_ ) );
            sres_   -= scalar * max( eval( lhs_ ), eval( rhs_ ) );
            osres_  -= scalar * max( eval( lhs_ ), eval( rhs_ ) );
            refres_ -= scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= scalar * max( eval( lhs_ ), eval( orhs_ ) );
            odres_  -= scalar * max( eval( lhs_ ), eval( orhs_ ) );
            sres_   -= scalar * max( eval( lhs_ ), eval( orhs_ ) );
            osres_  -= scalar * max( eval( lhs_ ), eval( orhs_ ) );
            refres_ -= scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= scalar * max( eval( olhs_ ), eval( rhs_ ) );
            odres_  -= scalar * max( eval( olhs_ ), eval( rhs_ ) );
            sres_   -= scalar * max( eval( olhs_ ), eval( rhs_ ) );
            osres_  -= scalar * max( eval( olhs_ ), eval( rhs_ ) );
            refres_ -= scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= scalar * max( eval( olhs_ ), eval( orhs_ ) );
            odres_  -= scalar * max( eval( olhs_ ), eval( orhs_ ) );
            sres_   -= scalar * max( eval( olhs_ ), eval( orhs_ ) );
            osres_  -= scalar * max( eval( olhs_ ), eval( orhs_ ) );
            refres_ -= scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled maximum with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled maximum with subtraction assignment with the given matrices
      {
         test_  = "Scaled maximum with subtraction assignment with the given matrices (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= max( lhs_, rhs_ ) * scalar;
            odres_  -= max( lhs_, rhs_ ) * scalar;
            sres_   -= max( lhs_, rhs_ ) * scalar;
            osres_  -= max( lhs_, rhs_ ) * scalar;
            refres_ -= ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= max( lhs_, orhs_ ) * scalar;
            odres_  -= max( lhs_, orhs_ ) * scalar;
            sres_   -= max( lhs_, orhs_ ) * scalar;
            osres_  -= max( lhs_, orhs_ ) * scalar;
            refres_ -= ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= max( olhs_, rhs_ ) * scalar;
            odres_  -= max( olhs_, rhs_ ) * scalar;
            sres_   -= max( olhs_, rhs_ ) * scalar;
            osres_  -= max( olhs_, rhs_ ) * scalar;
            refres_ -= ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= max( olhs_, orhs_ ) * scalar;
            odres_  -= max( olhs_, orhs_ ) * scalar;
            sres_   -= max( olhs_, orhs_ ) * scalar;
            osres_  -= max( olhs_, orhs_ ) * scalar;
            refres_ -= ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Scaled maximum with subtraction assignment with evaluated matrices (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            odres_  -= max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            sres_   -= max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            osres_  -= max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            refres_ -= eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            odres_  -= max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            sres_   -= max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            osres_  -= max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            refres_ -= eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            odres_  -= max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            sres_   -= max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            osres_  -= max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            refres_ -= eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            odres_  -= max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            sres_   -= max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            osres_  -= max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            refres_ -= eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled maximum with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled maximum with subtraction assignment with the given matrices
      {
         test_  = "Scaled maximum with subtraction assignment with the given matrices (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= max( lhs_, rhs_ ) / scalar;
            odres_  -= max( lhs_, rhs_ ) / scalar;
            sres_   -= max( lhs_, rhs_ ) / scalar;
            osres_  -= max( lhs_, rhs_ ) / scalar;
            refres_ -= ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= max( lhs_, orhs_ ) / scalar;
            odres_  -= max( lhs_, orhs_ ) / scalar;
            sres_   -= max( lhs_, orhs_ ) / scalar;
            osres_  -= max( lhs_, orhs_ ) / scalar;
            refres_ -= ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= max( olhs_, rhs_ ) / scalar;
            odres_  -= max( olhs_, rhs_ ) / scalar;
            sres_   -= max( olhs_, rhs_ ) / scalar;
            osres_  -= max( olhs_, rhs_ ) / scalar;
            refres_ -= ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= max( olhs_, orhs_ ) / scalar;
            odres_  -= max( olhs_, orhs_ ) / scalar;
            sres_   -= max( olhs_, orhs_ ) / scalar;
            osres_  -= max( olhs_, orhs_ ) / scalar;
            refres_ -= ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Scaled maximum with subtraction assignment with evaluated matrices (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            odres_  -= max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            sres_   -= max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            osres_  -= max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            refres_ -= eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            odres_  -= max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            sres_   -= max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            osres_  -= max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            refres_ -= eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            odres_  -= max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            sres_   -= max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            osres_  -= max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            refres_ -= eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            odres_  -= max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            sres_   -= max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            osres_  -= max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            refres_ -= eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled maximum with Schur product assignment (s*OP)
      //=====================================================================================

      // Scaled maximum with Schur product assignment with the given matrices
      {
         test_  = "Scaled maximum with Schur product assignment with the given matrices (s*OP)";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= scalar * max( lhs_, rhs_ );
            odres_  %= scalar * max( lhs_, rhs_ );
            sres_   %= scalar * max( lhs_, rhs_ );
            osres_  %= scalar * max( lhs_, rhs_ );
            refres_ %= scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= scalar * max( lhs_, orhs_ );
            odres_  %= scalar * max( lhs_, orhs_ );
            sres_   %= scalar * max( lhs_, orhs_ );
            osres_  %= scalar * max( lhs_, orhs_ );
            refres_ %= scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= scalar * max( olhs_, rhs_ );
            odres_  %= scalar * max( olhs_, rhs_ );
            sres_   %= scalar * max( olhs_, rhs_ );
            osres_  %= scalar * max( olhs_, rhs_ );
            refres_ %= scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= scalar * max( olhs_, orhs_ ) ;
            odres_  %= scalar * max( olhs_, orhs_ );
            sres_   %= scalar * max( olhs_, orhs_ );
            osres_  %= scalar * max( olhs_, orhs_ );
            refres_ %= scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled maximum with Schur product assignment with evaluated matrices
      {
         test_  = "Scaled maximum with Schur product assignment with evaluated matrices (s*OP)";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= scalar * max( eval( lhs_ ), eval( rhs_ ) );
            odres_  %= scalar * max( eval( lhs_ ), eval( rhs_ ) );
            sres_   %= scalar * max( eval( lhs_ ), eval( rhs_ ) );
            osres_  %= scalar * max( eval( lhs_ ), eval( rhs_ ) );
            refres_ %= scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= scalar * max( eval( lhs_ ), eval( orhs_ ) );
            odres_  %= scalar * max( eval( lhs_ ), eval( orhs_ ) );
            sres_   %= scalar * max( eval( lhs_ ), eval( orhs_ ) );
            osres_  %= scalar * max( eval( lhs_ ), eval( orhs_ ) );
            refres_ %= scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= scalar * max( eval( olhs_ ), eval( rhs_ ) );
            odres_  %= scalar * max( eval( olhs_ ), eval( rhs_ ) );
            sres_   %= scalar * max( eval( olhs_ ), eval( rhs_ ) );
            osres_  %= scalar * max( eval( olhs_ ), eval( rhs_ ) );
            refres_ %= scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= scalar * max( eval( olhs_ ), eval( orhs_ ) );
            odres_  %= scalar * max( eval( olhs_ ), eval( orhs_ ) );
            sres_   %= scalar * max( eval( olhs_ ), eval( orhs_ ) );
            osres_  %= scalar * max( eval( olhs_ ), eval( orhs_ ) );
            refres_ %= scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled maximum with Schur product assignment (OP*s)
      //=====================================================================================

      // Scaled maximum with Schur product assignment with the given matrices
      {
         test_  = "Scaled maximum with Schur product assignment with the given matrices (OP*s)";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= max( lhs_, rhs_ ) * scalar;
            odres_  %= max( lhs_, rhs_ ) * scalar;
            sres_   %= max( lhs_, rhs_ ) * scalar;
            osres_  %= max( lhs_, rhs_ ) * scalar;
            refres_ %= ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= max( lhs_, orhs_ ) * scalar;
            odres_  %= max( lhs_, orhs_ ) * scalar;
            sres_   %= max( lhs_, orhs_ ) * scalar;
            osres_  %= max( lhs_, orhs_ ) * scalar;
            refres_ %= ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= max( olhs_, rhs_ ) * scalar;
            odres_  %= max( olhs_, rhs_ ) * scalar;
            sres_   %= max( olhs_, rhs_ ) * scalar;
            osres_  %= max( olhs_, rhs_ ) * scalar;
            refres_ %= ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= max( olhs_, orhs_ ) * scalar;
            odres_  %= max( olhs_, orhs_ ) * scalar;
            sres_   %= max( olhs_, orhs_ ) * scalar;
            osres_  %= max( olhs_, orhs_ ) * scalar;
            refres_ %= ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled maximum with Schur product assignment with evaluated matrices
      {
         test_  = "Scaled maximum with Schur product assignment with evaluated matrices (OP*s)";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            odres_  %= max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            sres_   %= max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            osres_  %= max( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            refres_ %= eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            odres_  %= max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            sres_   %= max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            osres_  %= max( eval( lhs_ ), eval( orhs_ ) ) * scalar;
            refres_ %= eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            odres_  %= max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            sres_   %= max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            osres_  %= max( eval( olhs_ ), eval( rhs_ ) ) * scalar;
            refres_ %= eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            odres_  %= max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            sres_   %= max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            osres_  %= max( eval( olhs_ ), eval( orhs_ ) ) * scalar;
            refres_ %= eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled maximum with Schur product assignment (OP/s)
      //=====================================================================================

      // Scaled maximum with Schur product assignment with the given matrices
      {
         test_  = "Scaled maximum with Schur product assignment with the given matrices (OP/s)";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= max( lhs_, rhs_ ) / scalar;
            odres_  %= max( lhs_, rhs_ ) / scalar;
            sres_   %= max( lhs_, rhs_ ) / scalar;
            osres_  %= max( lhs_, rhs_ ) / scalar;
            refres_ %= ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= max( lhs_, orhs_ ) / scalar;
            odres_  %= max( lhs_, orhs_ ) / scalar;
            sres_   %= max( lhs_, orhs_ ) / scalar;
            osres_  %= max( lhs_, orhs_ ) / scalar;
            refres_ %= ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= max( olhs_, rhs_ ) / scalar;
            odres_  %= max( olhs_, rhs_ ) / scalar;
            sres_   %= max( olhs_, rhs_ ) / scalar;
            osres_  %= max( olhs_, rhs_ ) / scalar;
            refres_ %= ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= max( olhs_, orhs_ ) / scalar;
            odres_  %= max( olhs_, orhs_ ) / scalar;
            sres_   %= max( olhs_, orhs_ ) / scalar;
            osres_  %= max( olhs_, orhs_ ) / scalar;
            refres_ %= ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled maximum with Schur product assignment with evaluated matrices
      {
         test_  = "Scaled maximum with Schur product assignment with evaluated matrices (OP/s)";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            odres_  %= max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            sres_   %= max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            osres_  %= max( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            refres_ %= eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            odres_  %= max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            sres_   %= max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            osres_  %= max( eval( lhs_ ), eval( orhs_ ) ) / scalar;
            refres_ %= eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            odres_  %= max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            sres_   %= max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            osres_  %= max( eval( olhs_ ), eval( rhs_ ) ) / scalar;
            refres_ %= eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            odres_  %= max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            sres_   %= max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            osres_  %= max( eval( olhs_ ), eval( orhs_ ) ) / scalar;
            refres_ %= eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the transpose dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the transpose matrix maximum with plain assignment. In case any error
// resulting from the maximum or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANS_OPERATION > 1 )
   {
      //=====================================================================================
      // Transpose maximum
      //=====================================================================================

      // Transpose maximum with the given matrices
      {
         test_  = "Transpose maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initTransposeResults();
            tdres_  = trans( max( lhs_, rhs_ ) );
            todres_ = trans( max( lhs_, rhs_ ) );
            tsres_  = trans( max( lhs_, rhs_ ) );
            tosres_ = trans( max( lhs_, rhs_ ) );
            refres_ = trans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkTransposeResults<MT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = trans( max( lhs_, orhs_ ) );
            todres_ = trans( max( lhs_, orhs_ ) );
            tsres_  = trans( max( lhs_, orhs_ ) );
            tosres_ = trans( max( lhs_, orhs_ ) );
            refres_ = trans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkTransposeResults<MT1,OMT2>();

         try {
            initTransposeResults();
            tdres_  = trans( max( olhs_, rhs_ ) );
            todres_ = trans( max( olhs_, rhs_ ) );
            tsres_  = trans( max( olhs_, rhs_ ) );
            tosres_ = trans( max( olhs_, rhs_ ) );
            refres_ = trans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkTransposeResults<OMT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = trans( max( olhs_, orhs_ ) );
            todres_ = trans( max( olhs_, orhs_ ) );
            tsres_  = trans( max( olhs_, orhs_ ) );
            tosres_ = trans( max( olhs_, orhs_ ) );
            refres_ = trans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkTransposeResults<OMT1,OMT2>();
      }

      // Transpose maximum with evaluated matrices
      {
         test_  = "Transpose maximum with evaluated matrices";
         error_ = "Failed maximum operation";

         try {
            initTransposeResults();
            tdres_  = trans( max( eval( lhs_ ), eval( rhs_ ) ) );
            todres_ = trans( max( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_  = trans( max( eval( lhs_ ), eval( rhs_ ) ) );
            tosres_ = trans( max( eval( lhs_ ), eval( rhs_ ) ) );
            refres_ = trans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkTransposeResults<MT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = trans( max( eval( lhs_ ), eval( orhs_ ) ) );
            todres_ = trans( max( eval( lhs_ ), eval( orhs_ ) ) );
            tsres_  = trans( max( eval( lhs_ ), eval( orhs_ ) ) );
            tosres_ = trans( max( eval( lhs_ ), eval( orhs_ ) ) );
            refres_ = trans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkTransposeResults<MT1,OMT2>();

         try {
            initTransposeResults();
            tdres_  = trans( max( eval( olhs_ ), eval( rhs_ ) ) );
            todres_ = trans( max( eval( olhs_ ), eval( rhs_ ) ) );
            tsres_  = trans( max( eval( olhs_ ), eval( rhs_ ) ) );
            tosres_ = trans( max( eval( olhs_ ), eval( rhs_ ) ) );
            refres_ = trans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkTransposeResults<OMT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = trans( max( eval( olhs_ ), eval( orhs_ ) ) );
            todres_ = trans( max( eval( olhs_ ), eval( orhs_ ) ) );
            tsres_  = trans( max( eval( olhs_ ), eval( orhs_ ) ) );
            tosres_ = trans( max( eval( olhs_ ), eval( orhs_ ) ) );
            refres_ = trans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkTransposeResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the conjugate transpose dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the conjugate transpose matrix maximum with plain assignment. In case
// any error resulting from the maximum operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testCTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION > 1 )
   {
      //=====================================================================================
      // Conjugate transpose maximum
      //=====================================================================================

      // Conjugate transpose maximum with the given matrices
      {
         test_  = "Conjugate transpose maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initTransposeResults();
            tdres_  = ctrans( max( lhs_, rhs_ ) );
            todres_ = ctrans( max( lhs_, rhs_ ) );
            tsres_  = ctrans( max( lhs_, rhs_ ) );
            tosres_ = ctrans( max( lhs_, rhs_ ) );
            refres_ = ctrans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkTransposeResults<MT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = ctrans( max( lhs_, orhs_ ) );
            todres_ = ctrans( max( lhs_, orhs_ ) );
            tsres_  = ctrans( max( lhs_, orhs_ ) );
            tosres_ = ctrans( max( lhs_, orhs_ ) );
            refres_ = ctrans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkTransposeResults<MT1,OMT2>();

         try {
            initTransposeResults();
            tdres_  = ctrans( max( olhs_, rhs_ ) );
            todres_ = ctrans( max( olhs_, rhs_ ) );
            tsres_  = ctrans( max( olhs_, rhs_ ) );
            tosres_ = ctrans( max( olhs_, rhs_ ) );
            refres_ = ctrans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkTransposeResults<OMT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = ctrans( max( olhs_, orhs_ ) );
            todres_ = ctrans( max( olhs_, orhs_ ) );
            tsres_  = ctrans( max( olhs_, orhs_ ) );
            tosres_ = ctrans( max( olhs_, orhs_ ) );
            refres_ = ctrans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkTransposeResults<OMT1,OMT2>();
      }

      // Conjugate transpose maximum with evaluated matrices
      {
         test_  = "Conjugate transpose maximum with evaluated matrices";
         error_ = "Failed maximum operation";

         try {
            initTransposeResults();
            tdres_  = ctrans( max( eval( lhs_ ), eval( rhs_ ) ) );
            todres_ = ctrans( max( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_  = ctrans( max( eval( lhs_ ), eval( rhs_ ) ) );
            tosres_ = ctrans( max( eval( lhs_ ), eval( rhs_ ) ) );
            refres_ = ctrans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkTransposeResults<MT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = ctrans( max( eval( lhs_ ), eval( orhs_ ) ) );
            todres_ = ctrans( max( eval( lhs_ ), eval( orhs_ ) ) );
            tsres_  = ctrans( max( eval( lhs_ ), eval( orhs_ ) ) );
            tosres_ = ctrans( max( eval( lhs_ ), eval( orhs_ ) ) );
            refres_ = ctrans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkTransposeResults<MT1,OMT2>();

         try {
            initTransposeResults();
            tdres_  = ctrans( max( eval( olhs_ ), eval( rhs_ ) ) );
            todres_ = ctrans( max( eval( olhs_ ), eval( rhs_ ) ) );
            tsres_  = ctrans( max( eval( olhs_ ), eval( rhs_ ) ) );
            tosres_ = ctrans( max( eval( olhs_ ), eval( rhs_ ) ) );
            refres_ = ctrans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkTransposeResults<OMT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = ctrans( max( eval( olhs_ ), eval( orhs_ ) ) );
            todres_ = ctrans( max( eval( olhs_ ), eval( orhs_ ) ) );
            tsres_  = ctrans( max( eval( olhs_ ), eval( orhs_ ) ) );
            tosres_ = ctrans( max( eval( olhs_ ), eval( orhs_ ) ) );
            refres_ = ctrans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkTransposeResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the abs dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the abs matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testAbsOperation()
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
/*!\brief Testing the conjugate dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the conjugate matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testConjOperation()
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
/*!\brief Testing the \a real dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the \a real matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testRealOperation()
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
/*!\brief Testing the \a imag dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the \a imag matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testImagOperation()
{
#if BLAZETEST_MATHTEST_TEST_IMAG_OPERATION
   if( BLAZETEST_MATHTEST_TEST_IMAG_OPERATION > 1 &&
       ( !blaze::IsHermitian<DRE>::value || blaze::isSymmetric( imag( max( lhs_, rhs_ ) ) ) ) )
   {
      testCustomOperation( blaze::Imag(), "imag" );
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the \a inv dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the \a inv matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testInvOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_INV_OPERATION && BLAZETEST_MATHTEST_LAPACK_MODE
   if( BLAZETEST_MATHTEST_TEST_INV_OPERATION > 1 )
   {
      if( !isSquare( max( lhs_, rhs_ ) ) || blaze::isDefault( det( max( lhs_, rhs_ ) ) ) )
         return;

      testCustomOperation( blaze::Inv(), "inv" );
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the \a inv dense matrix/dense matrix maximum.
//
// \return void
//
// This function is called in case the \a inv matrix/matrix maximum operation is not available
// for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testInvOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the evaluated dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the evaluated matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testEvalOperation()
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
/*!\brief Testing the serialized dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the serialized matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testSerialOperation()
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
/*!\brief Testing the non-aliased dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the non-aliased matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testNoAliasOperation()
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
/*!\brief Testing the non-SIMD dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the non-SIMD matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testNoSIMDOperation()
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
/*!\brief Testing the symmetric dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the symmetric matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testDeclSymOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION > 1 )
   {
      if( ( !blaze::IsDiagonal<MT1>::value && blaze::IsTriangular<MT1>::value ) ||
          ( !blaze::IsDiagonal<MT2>::value && blaze::IsTriangular<MT2>::value ) ||
          ( !blaze::IsDiagonal<MT1>::value && blaze::IsHermitian<MT1>::value && blaze::IsComplex<ET1>::value ) ||
          ( !blaze::IsDiagonal<MT2>::value && blaze::IsHermitian<MT2>::value && blaze::IsComplex<ET2>::value ) ||
          ( lhs_.rows() != lhs_.columns() ) )
         return;


      //=====================================================================================
      // Test-specific setup of the left-hand side operand
      //=====================================================================================

      MT1  lhs ( lhs_ * trans( lhs_ ) );
      OMT1 olhs( lhs );


      //=====================================================================================
      // Test-specific setup of the right-hand side operand
      //=====================================================================================

      MT2  rhs ( rhs_ * trans( rhs_ ) );
      OMT2 orhs( rhs );


      //=====================================================================================
      // Test-specific setup of the reference matrix
      //=====================================================================================

      RT ref( max( lhs, rhs ) );


      //=====================================================================================
      // Declsym maximum
      //=====================================================================================

      // Declsym maximum with the given matrices
      {
         test_  = "Declsym maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = declsym( max( lhs, rhs ) );
            odres_  = declsym( max( lhs, rhs ) );
            sres_   = declsym( max( lhs, rhs ) );
            osres_  = declsym( max( lhs, rhs ) );
            refres_ = declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declsym( max( lhs, orhs ) );
            odres_  = declsym( max( lhs, orhs ) );
            sres_   = declsym( max( lhs, orhs ) );
            osres_  = declsym( max( lhs, orhs ) );
            refres_ = declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declsym( max( olhs, rhs ) );
            odres_  = declsym( max( olhs, rhs ) );
            sres_   = declsym( max( olhs, rhs ) );
            osres_  = declsym( max( olhs, rhs ) );
            refres_ = declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declsym( max( olhs, orhs ) );
            odres_  = declsym( max( olhs, orhs ) );
            sres_   = declsym( max( olhs, orhs ) );
            osres_  = declsym( max( olhs, orhs ) );
            refres_ = declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declsym maximum with evaluated matrices
      {
         test_  = "Declsym maximum with evaluated left-hand side matrix";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = declsym( max( eval( lhs ), eval( rhs ) ) );
            odres_  = declsym( max( eval( lhs ), eval( rhs ) ) );
            sres_   = declsym( max( eval( lhs ), eval( rhs ) ) );
            osres_  = declsym( max( eval( lhs ), eval( rhs ) ) );
            refres_ = declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declsym( max( eval( lhs ), eval( orhs ) ) );
            odres_  = declsym( max( eval( lhs ), eval( orhs ) ) );
            sres_   = declsym( max( eval( lhs ), eval( orhs ) ) );
            osres_  = declsym( max( eval( lhs ), eval( orhs ) ) );
            refres_ = declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declsym( max( eval( olhs ), eval( rhs ) ) );
            odres_  = declsym( max( eval( olhs ), eval( rhs ) ) );
            sres_   = declsym( max( eval( olhs ), eval( rhs ) ) );
            osres_  = declsym( max( eval( olhs ), eval( rhs ) ) );
            refres_ = declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declsym( max( eval( olhs ), eval( orhs ) ) );
            odres_  = declsym( max( eval( olhs ), eval( orhs ) ) );
            sres_   = declsym( max( eval( olhs ), eval( orhs ) ) );
            osres_  = declsym( max( eval( olhs ), eval( orhs ) ) );
            refres_ = declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declsym maximum with addition assignment
      //=====================================================================================

      // Declsym maximum with addition assignment with the given matrices
      {
         test_  = "Declsym maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declsym( max( lhs, rhs ) );
            odres_  += declsym( max( lhs, rhs ) );
            sres_   += declsym( max( lhs, rhs ) );
            osres_  += declsym( max( lhs, rhs ) );
            refres_ += declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declsym( max( lhs, orhs ) );
            odres_  += declsym( max( lhs, orhs ) );
            sres_   += declsym( max( lhs, orhs ) );
            osres_  += declsym( max( lhs, orhs ) );
            refres_ += declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declsym( max( olhs, rhs ) );
            odres_  += declsym( max( olhs, rhs ) );
            sres_   += declsym( max( olhs, rhs ) );
            osres_  += declsym( max( olhs, rhs ) );
            refres_ += declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declsym( max( olhs, orhs ) );
            odres_  += declsym( max( olhs, orhs ) );
            sres_   += declsym( max( olhs, orhs ) );
            osres_  += declsym( max( olhs, orhs ) );
            refres_ += declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declsym maximum with addition assignment with evaluated matrices
      {
         test_  = "Declsym maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declsym( max( eval( lhs ), eval( rhs ) ) );
            odres_  += declsym( max( eval( lhs ), eval( rhs ) ) );
            sres_   += declsym( max( eval( lhs ), eval( rhs ) ) );
            osres_  += declsym( max( eval( lhs ), eval( rhs ) ) );
            refres_ += declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declsym( max( eval( lhs ), eval( orhs ) ) );
            odres_  += declsym( max( eval( lhs ), eval( orhs ) ) );
            sres_   += declsym( max( eval( lhs ), eval( orhs ) ) );
            osres_  += declsym( max( eval( lhs ), eval( orhs ) ) );
            refres_ += declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declsym( max( eval( olhs ), eval( rhs ) ) );
            odres_  += declsym( max( eval( olhs ), eval( rhs ) ) );
            sres_   += declsym( max( eval( olhs ), eval( rhs ) ) );
            osres_  += declsym( max( eval( olhs ), eval( rhs ) ) );
            refres_ += declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declsym( max( eval( olhs ), eval( orhs ) ) );
            odres_  += declsym( max( eval( olhs ), eval( orhs ) ) );
            sres_   += declsym( max( eval( olhs ), eval( orhs ) ) );
            osres_  += declsym( max( eval( olhs ), eval( orhs ) ) );
            refres_ += declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declsym maximum with subtraction assignment
      //=====================================================================================

      // Declsym maximum with subtraction assignment with the given matrices
      {
         test_  = "Declsym maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declsym( max( lhs, rhs ) );
            odres_  -= declsym( max( lhs, rhs ) );
            sres_   -= declsym( max( lhs, rhs ) );
            osres_  -= declsym( max( lhs, rhs ) );
            refres_ -= declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declsym( max( lhs, orhs ) );
            odres_  -= declsym( max( lhs, orhs ) );
            sres_   -= declsym( max( lhs, orhs ) );
            osres_  -= declsym( max( lhs, orhs ) );
            refres_ -= declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declsym( max( olhs, rhs ) );
            odres_  -= declsym( max( olhs, rhs ) );
            sres_   -= declsym( max( olhs, rhs ) );
            osres_  -= declsym( max( olhs, rhs ) );
            refres_ -= declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declsym( max( olhs, orhs ) );
            odres_  -= declsym( max( olhs, orhs ) );
            sres_   -= declsym( max( olhs, orhs ) );
            osres_  -= declsym( max( olhs, orhs ) );
            refres_ -= declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declsym maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Declsym maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declsym( max( eval( lhs ), eval( rhs ) ) );
            odres_  -= declsym( max( eval( lhs ), eval( rhs ) ) );
            sres_   -= declsym( max( eval( lhs ), eval( rhs ) ) );
            osres_  -= declsym( max( eval( lhs ), eval( rhs ) ) );
            refres_ -= declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declsym( max( eval( lhs ), eval( orhs ) ) );
            odres_  -= declsym( max( eval( lhs ), eval( orhs ) ) );
            sres_   -= declsym( max( eval( lhs ), eval( orhs ) ) );
            osres_  -= declsym( max( eval( lhs ), eval( orhs ) ) );
            refres_ -= declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declsym( max( eval( olhs ), eval( rhs ) ) );
            odres_  -= declsym( max( eval( olhs ), eval( rhs ) ) );
            sres_   -= declsym( max( eval( olhs ), eval( rhs ) ) );
            osres_  -= declsym( max( eval( olhs ), eval( rhs ) ) );
            refres_ -= declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declsym( max( eval( olhs ), eval( orhs ) ) );
            odres_  -= declsym( max( eval( olhs ), eval( orhs ) ) );
            sres_   -= declsym( max( eval( olhs ), eval( orhs ) ) );
            osres_  -= declsym( max( eval( olhs ), eval( orhs ) ) );
            refres_ -= declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declsym maximum with Schur product assignment
      //=====================================================================================

      // Declsym maximum with Schur product assignment with the given matrices
      {
         test_  = "Declsym maximum with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declsym( max( lhs, rhs ) );
            odres_  %= declsym( max( lhs, rhs ) );
            sres_   %= declsym( max( lhs, rhs ) );
            osres_  %= declsym( max( lhs, rhs ) );
            refres_ %= declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declsym( max( lhs, orhs ) );
            odres_  %= declsym( max( lhs, orhs ) );
            sres_   %= declsym( max( lhs, orhs ) );
            osres_  %= declsym( max( lhs, orhs ) );
            refres_ %= declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declsym( max( olhs, rhs ) );
            odres_  %= declsym( max( olhs, rhs ) );
            sres_   %= declsym( max( olhs, rhs ) );
            osres_  %= declsym( max( olhs, rhs ) );
            refres_ %= declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declsym( max( olhs, orhs ) );
            odres_  %= declsym( max( olhs, orhs ) );
            sres_   %= declsym( max( olhs, orhs ) );
            osres_  %= declsym( max( olhs, orhs ) );
            refres_ %= declsym( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declsym maximum with Schur product assignment with evaluated matrices
      {
         test_  = "Declsym maximum with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declsym( max( eval( lhs ), eval( rhs ) ) );
            odres_  %= declsym( max( eval( lhs ), eval( rhs ) ) );
            sres_   %= declsym( max( eval( lhs ), eval( rhs ) ) );
            osres_  %= declsym( max( eval( lhs ), eval( rhs ) ) );
            refres_ %= declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declsym( max( eval( lhs ), eval( orhs ) ) );
            odres_  %= declsym( max( eval( lhs ), eval( orhs ) ) );
            sres_   %= declsym( max( eval( lhs ), eval( orhs ) ) );
            osres_  %= declsym( max( eval( lhs ), eval( orhs ) ) );
            refres_ %= declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declsym( max( eval( olhs ), eval( rhs ) ) );
            odres_  %= declsym( max( eval( olhs ), eval( rhs ) ) );
            sres_   %= declsym( max( eval( olhs ), eval( rhs ) ) );
            osres_  %= declsym( max( eval( olhs ), eval( rhs ) ) );
            refres_ %= declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declsym( max( eval( olhs ), eval( orhs ) ) );
            odres_  %= declsym( max( eval( olhs ), eval( orhs ) ) );
            sres_   %= declsym( max( eval( olhs ), eval( orhs ) ) );
            osres_  %= declsym( max( eval( olhs ), eval( orhs ) ) );
            refres_ %= declsym( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the symmetric dense matrix/dense matrix maximum operation.
//
// \return void
//
// This function is called in case the symmetric matrix/matrix maximum operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testDeclSymOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the Hermitian dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the Hermitian matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testDeclHermOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION > 1 )
   {
      if( ( !blaze::IsDiagonal<MT1>::value && blaze::IsTriangular<MT1>::value ) ||
          ( !blaze::IsDiagonal<MT2>::value && blaze::IsTriangular<MT2>::value ) ||
          ( !blaze::IsDiagonal<MT1>::value && blaze::IsSymmetric<MT1>::value && blaze::IsComplex<ET1>::value ) ||
          ( !blaze::IsDiagonal<MT2>::value && blaze::IsSymmetric<MT2>::value && blaze::IsComplex<ET2>::value ) ||
          ( lhs_.rows() != lhs_.columns() ) )
         return;


      //=====================================================================================
      // Test-specific setup of the left-hand side operand
      //=====================================================================================

      MT1  lhs ( lhs_ * ctrans( lhs_ ) );
      OMT1 olhs( lhs );


      //=====================================================================================
      // Test-specific setup of the right-hand side operand
      //=====================================================================================

      MT2  rhs ( rhs_ * ctrans( rhs_ ) );
      OMT2 orhs( rhs );


      //=====================================================================================
      // Test-specific setup of the reference matrix
      //=====================================================================================

      RT ref( max( lhs, rhs ) );


      //=====================================================================================
      // Declherm maximum
      //=====================================================================================

      // Declherm maximum with the given matrices
      {
         test_  = "Declherm maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = declherm( max( lhs, rhs ) );
            odres_  = declherm( max( lhs, rhs ) );
            sres_   = declherm( max( lhs, rhs ) );
            osres_  = declherm( max( lhs, rhs ) );
            refres_ = declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declherm( max( lhs, orhs ) );
            odres_  = declherm( max( lhs, orhs ) );
            sres_   = declherm( max( lhs, orhs ) );
            osres_  = declherm( max( lhs, orhs ) );
            refres_ = declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declherm( max( olhs, rhs ) );
            odres_  = declherm( max( olhs, rhs ) );
            sres_   = declherm( max( olhs, rhs ) );
            osres_  = declherm( max( olhs, rhs ) );
            refres_ = declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declherm( max( olhs, orhs ) );
            odres_  = declherm( max( olhs, orhs ) );
            sres_   = declherm( max( olhs, orhs ) );
            osres_  = declherm( max( olhs, orhs ) );
            refres_ = declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declherm maximum with evaluated matrices
      {
         test_  = "Declherm maximum with evaluated left-hand side matrix";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = declherm( max( eval( lhs ), eval( rhs ) ) );
            odres_  = declherm( max( eval( lhs ), eval( rhs ) ) );
            sres_   = declherm( max( eval( lhs ), eval( rhs ) ) );
            osres_  = declherm( max( eval( lhs ), eval( rhs ) ) );
            refres_ = declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declherm( max( eval( lhs ), eval( orhs ) ) );
            odres_  = declherm( max( eval( lhs ), eval( orhs ) ) );
            sres_   = declherm( max( eval( lhs ), eval( orhs ) ) );
            osres_  = declherm( max( eval( lhs ), eval( orhs ) ) );
            refres_ = declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declherm( max( eval( olhs ), eval( rhs ) ) );
            odres_  = declherm( max( eval( olhs ), eval( rhs ) ) );
            sres_   = declherm( max( eval( olhs ), eval( rhs ) ) );
            osres_  = declherm( max( eval( olhs ), eval( rhs ) ) );
            refres_ = declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declherm( max( eval( olhs ), eval( orhs ) ) );
            odres_  = declherm( max( eval( olhs ), eval( orhs ) ) );
            sres_   = declherm( max( eval( olhs ), eval( orhs ) ) );
            osres_  = declherm( max( eval( olhs ), eval( orhs ) ) );
            refres_ = declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declherm maximum with addition assignment
      //=====================================================================================

      // Declherm maximum with addition assignment with the given matrices
      {
         test_  = "Declherm maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declherm( max( lhs, rhs ) );
            odres_  += declherm( max( lhs, rhs ) );
            sres_   += declherm( max( lhs, rhs ) );
            osres_  += declherm( max( lhs, rhs ) );
            refres_ += declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declherm( max( lhs, orhs ) );
            odres_  += declherm( max( lhs, orhs ) );
            sres_   += declherm( max( lhs, orhs ) );
            osres_  += declherm( max( lhs, orhs ) );
            refres_ += declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declherm( max( olhs, rhs ) );
            odres_  += declherm( max( olhs, rhs ) );
            sres_   += declherm( max( olhs, rhs ) );
            osres_  += declherm( max( olhs, rhs ) );
            refres_ += declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declherm( max( olhs, orhs ) );
            odres_  += declherm( max( olhs, orhs ) );
            sres_   += declherm( max( olhs, orhs ) );
            osres_  += declherm( max( olhs, orhs ) );
            refres_ += declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declherm maximum with addition assignment with evaluated matrices
      {
         test_  = "Declherm maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declherm( max( eval( lhs ), eval( rhs ) ) );
            odres_  += declherm( max( eval( lhs ), eval( rhs ) ) );
            sres_   += declherm( max( eval( lhs ), eval( rhs ) ) );
            osres_  += declherm( max( eval( lhs ), eval( rhs ) ) );
            refres_ += declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declherm( max( eval( lhs ), eval( orhs ) ) );
            odres_  += declherm( max( eval( lhs ), eval( orhs ) ) );
            sres_   += declherm( max( eval( lhs ), eval( orhs ) ) );
            osres_  += declherm( max( eval( lhs ), eval( orhs ) ) );
            refres_ += declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declherm( max( eval( olhs ), eval( rhs ) ) );
            odres_  += declherm( max( eval( olhs ), eval( rhs ) ) );
            sres_   += declherm( max( eval( olhs ), eval( rhs ) ) );
            osres_  += declherm( max( eval( olhs ), eval( rhs ) ) );
            refres_ += declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declherm( max( eval( olhs ), eval( orhs ) ) );
            odres_  += declherm( max( eval( olhs ), eval( orhs ) ) );
            sres_   += declherm( max( eval( olhs ), eval( orhs ) ) );
            osres_  += declherm( max( eval( olhs ), eval( orhs ) ) );
            refres_ += declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declherm maximum with subtraction assignment
      //=====================================================================================

      // Declherm maximum with subtraction assignment with the given matrices
      {
         test_  = "Declherm maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declherm( max( lhs, rhs ) );
            odres_  -= declherm( max( lhs, rhs ) );
            sres_   -= declherm( max( lhs, rhs ) );
            osres_  -= declherm( max( lhs, rhs ) );
            refres_ -= declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declherm( max( lhs, orhs ) );
            odres_  -= declherm( max( lhs, orhs ) );
            sres_   -= declherm( max( lhs, orhs ) );
            osres_  -= declherm( max( lhs, orhs ) );
            refres_ -= declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declherm( max( olhs, rhs ) );
            odres_  -= declherm( max( olhs, rhs ) );
            sres_   -= declherm( max( olhs, rhs ) );
            osres_  -= declherm( max( olhs, rhs ) );
            refres_ -= declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declherm( max( olhs, orhs ) );
            odres_  -= declherm( max( olhs, orhs ) );
            sres_   -= declherm( max( olhs, orhs ) );
            osres_  -= declherm( max( olhs, orhs ) );
            refres_ -= declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declherm maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Declherm maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declherm( max( eval( lhs ), eval( rhs ) ) );
            odres_  -= declherm( max( eval( lhs ), eval( rhs ) ) );
            sres_   -= declherm( max( eval( lhs ), eval( rhs ) ) );
            osres_  -= declherm( max( eval( lhs ), eval( rhs ) ) );
            refres_ -= declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declherm( max( eval( lhs ), eval( orhs ) ) );
            odres_  -= declherm( max( eval( lhs ), eval( orhs ) ) );
            sres_   -= declherm( max( eval( lhs ), eval( orhs ) ) );
            osres_  -= declherm( max( eval( lhs ), eval( orhs ) ) );
            refres_ -= declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declherm( max( eval( olhs ), eval( rhs ) ) );
            odres_  -= declherm( max( eval( olhs ), eval( rhs ) ) );
            sres_   -= declherm( max( eval( olhs ), eval( rhs ) ) );
            osres_  -= declherm( max( eval( olhs ), eval( rhs ) ) );
            refres_ -= declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declherm( max( eval( olhs ), eval( orhs ) ) );
            odres_  -= declherm( max( eval( olhs ), eval( orhs ) ) );
            sres_   -= declherm( max( eval( olhs ), eval( orhs ) ) );
            osres_  -= declherm( max( eval( olhs ), eval( orhs ) ) );
            refres_ -= declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declherm maximum with Schur product assignment
      //=====================================================================================

      // Declherm maximum with Schur product assignment with the given matrices
      {
         test_  = "Declherm maximum with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declherm( max( lhs, rhs ) );
            odres_  %= declherm( max( lhs, rhs ) );
            sres_   %= declherm( max( lhs, rhs ) );
            osres_  %= declherm( max( lhs, rhs ) );
            refres_ %= declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declherm( max( lhs, orhs ) );
            odres_  %= declherm( max( lhs, orhs ) );
            sres_   %= declherm( max( lhs, orhs ) );
            osres_  %= declherm( max( lhs, orhs ) );
            refres_ %= declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declherm( max( olhs, rhs ) );
            odres_  %= declherm( max( olhs, rhs ) );
            sres_   %= declherm( max( olhs, rhs ) );
            osres_  %= declherm( max( olhs, rhs ) );
            refres_ %= declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declherm( max( olhs, orhs ) );
            odres_  %= declherm( max( olhs, orhs ) );
            sres_   %= declherm( max( olhs, orhs ) );
            osres_  %= declherm( max( olhs, orhs ) );
            refres_ %= declherm( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declherm maximum with Schur product assignment with evaluated matrices
      {
         test_  = "Declherm maximum with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declherm( max( eval( lhs ), eval( rhs ) ) );
            odres_  %= declherm( max( eval( lhs ), eval( rhs ) ) );
            sres_   %= declherm( max( eval( lhs ), eval( rhs ) ) );
            osres_  %= declherm( max( eval( lhs ), eval( rhs ) ) );
            refres_ %= declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declherm( max( eval( lhs ), eval( orhs ) ) );
            odres_  %= declherm( max( eval( lhs ), eval( orhs ) ) );
            sres_   %= declherm( max( eval( lhs ), eval( orhs ) ) );
            osres_  %= declherm( max( eval( lhs ), eval( orhs ) ) );
            refres_ %= declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declherm( max( eval( olhs ), eval( rhs ) ) );
            odres_  %= declherm( max( eval( olhs ), eval( rhs ) ) );
            sres_   %= declherm( max( eval( olhs ), eval( rhs ) ) );
            osres_  %= declherm( max( eval( olhs ), eval( rhs ) ) );
            refres_ %= declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declherm( max( eval( olhs ), eval( orhs ) ) );
            odres_  %= declherm( max( eval( olhs ), eval( orhs ) ) );
            sres_   %= declherm( max( eval( olhs ), eval( orhs ) ) );
            osres_  %= declherm( max( eval( olhs ), eval( orhs ) ) );
            refres_ %= declherm( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the Hermitian dense matrix/dense matrix maximum operation.
//
// \return void
//
// This function is called in case the Hermitian matrix/matrix maximum operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testDeclHermOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the lower dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the lower matrix maximum with plain assignment, addition assignment,
// and subtraction assignment. In case any error resulting from the maximum operation or the
// subsequent assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testDeclLowOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLLOW_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLLOW_OPERATION > 1 )
   {
      if( lhs_.rows() != lhs_.columns() )
         return;


      //=====================================================================================
      // Test-specific setup of the left-hand side operand
      //=====================================================================================

      MT1 lhs( lhs_ );

      blaze::resetUpper( lhs );

      OMT1 olhs( lhs );


      //=====================================================================================
      // Test-specific setup of the right-hand side operand
      //=====================================================================================

      MT2 rhs( rhs_ );

      blaze::resetUpper( rhs );

      OMT2 orhs( rhs );


      //=====================================================================================
      // Test-specific setup of the reference matrix
      //=====================================================================================

      RT ref( max( lhs, rhs ) );


      //=====================================================================================
      // Decllow maximum
      //=====================================================================================

      // Decllow maximum with the given matrices
      {
         test_  = "Decllow maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = decllow( max( lhs, rhs ) );
            odres_  = decllow( max( lhs, rhs ) );
            sres_   = decllow( max( lhs, rhs ) );
            osres_  = decllow( max( lhs, rhs ) );
            refres_ = decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = decllow( max( lhs, orhs ) );
            odres_  = decllow( max( lhs, orhs ) );
            sres_   = decllow( max( lhs, orhs ) );
            osres_  = decllow( max( lhs, orhs ) );
            refres_ = decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = decllow( max( olhs, rhs ) );
            odres_  = decllow( max( olhs, rhs ) );
            sres_   = decllow( max( olhs, rhs ) );
            osres_  = decllow( max( olhs, rhs ) );
            refres_ = decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = decllow( max( olhs, orhs ) );
            odres_  = decllow( max( olhs, orhs ) );
            sres_   = decllow( max( olhs, orhs ) );
            osres_  = decllow( max( olhs, orhs ) );
            refres_ = decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decllow maximum with evaluated matrices
      {
         test_  = "Decllow maximum with evaluated left-hand side matrix";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = decllow( max( eval( lhs ), eval( rhs ) ) );
            odres_  = decllow( max( eval( lhs ), eval( rhs ) ) );
            sres_   = decllow( max( eval( lhs ), eval( rhs ) ) );
            osres_  = decllow( max( eval( lhs ), eval( rhs ) ) );
            refres_ = decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = decllow( max( eval( lhs ), eval( orhs ) ) );
            odres_  = decllow( max( eval( lhs ), eval( orhs ) ) );
            sres_   = decllow( max( eval( lhs ), eval( orhs ) ) );
            osres_  = decllow( max( eval( lhs ), eval( orhs ) ) );
            refres_ = decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = decllow( max( eval( olhs ), eval( rhs ) ) );
            odres_  = decllow( max( eval( olhs ), eval( rhs ) ) );
            sres_   = decllow( max( eval( olhs ), eval( rhs ) ) );
            osres_  = decllow( max( eval( olhs ), eval( rhs ) ) );
            refres_ = decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = decllow( max( eval( olhs ), eval( orhs ) ) );
            odres_  = decllow( max( eval( olhs ), eval( orhs ) ) );
            sres_   = decllow( max( eval( olhs ), eval( orhs ) ) );
            osres_  = decllow( max( eval( olhs ), eval( orhs ) ) );
            refres_ = decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Decllow maximum with addition assignment
      //=====================================================================================

      // Decllow maximum with addition assignment with the given matrices
      {
         test_  = "Decllow maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += decllow( max( lhs, rhs ) );
            odres_  += decllow( max( lhs, rhs ) );
            sres_   += decllow( max( lhs, rhs ) );
            osres_  += decllow( max( lhs, rhs ) );
            refres_ += decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += decllow( max( lhs, orhs ) );
            odres_  += decllow( max( lhs, orhs ) );
            sres_   += decllow( max( lhs, orhs ) );
            osres_  += decllow( max( lhs, orhs ) );
            refres_ += decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += decllow( max( olhs, rhs ) );
            odres_  += decllow( max( olhs, rhs ) );
            sres_   += decllow( max( olhs, rhs ) );
            osres_  += decllow( max( olhs, rhs ) );
            refres_ += decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += decllow( max( olhs, orhs ) );
            odres_  += decllow( max( olhs, orhs ) );
            sres_   += decllow( max( olhs, orhs ) );
            osres_  += decllow( max( olhs, orhs ) );
            refres_ += decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decllow maximum with addition assignment with evaluated matrices
      {
         test_  = "Decllow maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += decllow( max( eval( lhs ), eval( rhs ) ) );
            odres_  += decllow( max( eval( lhs ), eval( rhs ) ) );
            sres_   += decllow( max( eval( lhs ), eval( rhs ) ) );
            osres_  += decllow( max( eval( lhs ), eval( rhs ) ) );
            refres_ += decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += decllow( max( eval( lhs ), eval( orhs ) ) );
            odres_  += decllow( max( eval( lhs ), eval( orhs ) ) );
            sres_   += decllow( max( eval( lhs ), eval( orhs ) ) );
            osres_  += decllow( max( eval( lhs ), eval( orhs ) ) );
            refres_ += decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += decllow( max( eval( olhs ), eval( rhs ) ) );
            odres_  += decllow( max( eval( olhs ), eval( rhs ) ) );
            sres_   += decllow( max( eval( olhs ), eval( rhs ) ) );
            osres_  += decllow( max( eval( olhs ), eval( rhs ) ) );
            refres_ += decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += decllow( max( eval( olhs ), eval( orhs ) ) );
            odres_  += decllow( max( eval( olhs ), eval( orhs ) ) );
            sres_   += decllow( max( eval( olhs ), eval( orhs ) ) );
            osres_  += decllow( max( eval( olhs ), eval( orhs ) ) );
            refres_ += decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Decllow maximum with subtraction assignment
      //=====================================================================================

      // Decllow maximum with subtraction assignment with the given matrices
      {
         test_  = "Decllow maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= decllow( max( lhs, rhs ) );
            odres_  -= decllow( max( lhs, rhs ) );
            sres_   -= decllow( max( lhs, rhs ) );
            osres_  -= decllow( max( lhs, rhs ) );
            refres_ -= decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= decllow( max( lhs, orhs ) );
            odres_  -= decllow( max( lhs, orhs ) );
            sres_   -= decllow( max( lhs, orhs ) );
            osres_  -= decllow( max( lhs, orhs ) );
            refres_ -= decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= decllow( max( olhs, rhs ) );
            odres_  -= decllow( max( olhs, rhs ) );
            sres_   -= decllow( max( olhs, rhs ) );
            osres_  -= decllow( max( olhs, rhs ) );
            refres_ -= decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= decllow( max( olhs, orhs ) );
            odres_  -= decllow( max( olhs, orhs ) );
            sres_   -= decllow( max( olhs, orhs ) );
            osres_  -= decllow( max( olhs, orhs ) );
            refres_ -= decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decllow maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Decllow maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= decllow( max( eval( lhs ), eval( rhs ) ) );
            odres_  -= decllow( max( eval( lhs ), eval( rhs ) ) );
            sres_   -= decllow( max( eval( lhs ), eval( rhs ) ) );
            osres_  -= decllow( max( eval( lhs ), eval( rhs ) ) );
            refres_ -= decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= decllow( max( eval( lhs ), eval( orhs ) ) );
            odres_  -= decllow( max( eval( lhs ), eval( orhs ) ) );
            sres_   -= decllow( max( eval( lhs ), eval( orhs ) ) );
            osres_  -= decllow( max( eval( lhs ), eval( orhs ) ) );
            refres_ -= decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= decllow( max( eval( olhs ), eval( rhs ) ) );
            odres_  -= decllow( max( eval( olhs ), eval( rhs ) ) );
            sres_   -= decllow( max( eval( olhs ), eval( rhs ) ) );
            osres_  -= decllow( max( eval( olhs ), eval( rhs ) ) );
            refres_ -= decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= decllow( max( eval( olhs ), eval( orhs ) ) );
            odres_  -= decllow( max( eval( olhs ), eval( orhs ) ) );
            sres_   -= decllow( max( eval( olhs ), eval( orhs ) ) );
            osres_  -= decllow( max( eval( olhs ), eval( orhs ) ) );
            refres_ -= decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Decllow maximum with Schur product assignment
      //=====================================================================================

      // Decllow maximum with Schur product assignment with the given matrices
      {
         test_  = "Decllow maximum with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= decllow( max( lhs, rhs ) );
            odres_  %= decllow( max( lhs, rhs ) );
            sres_   %= decllow( max( lhs, rhs ) );
            osres_  %= decllow( max( lhs, rhs ) );
            refres_ %= decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= decllow( max( lhs, orhs ) );
            odres_  %= decllow( max( lhs, orhs ) );
            sres_   %= decllow( max( lhs, orhs ) );
            osres_  %= decllow( max( lhs, orhs ) );
            refres_ %= decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= decllow( max( olhs, rhs ) );
            odres_  %= decllow( max( olhs, rhs ) );
            sres_   %= decllow( max( olhs, rhs ) );
            osres_  %= decllow( max( olhs, rhs ) );
            refres_ %= decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= decllow( max( olhs, orhs ) );
            odres_  %= decllow( max( olhs, orhs ) );
            sres_   %= decllow( max( olhs, orhs ) );
            osres_  %= decllow( max( olhs, orhs ) );
            refres_ %= decllow( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decllow maximum with Schur product assignment with evaluated matrices
      {
         test_  = "Decllow maximum with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= decllow( max( eval( lhs ), eval( rhs ) ) );
            odres_  %= decllow( max( eval( lhs ), eval( rhs ) ) );
            sres_   %= decllow( max( eval( lhs ), eval( rhs ) ) );
            osres_  %= decllow( max( eval( lhs ), eval( rhs ) ) );
            refres_ %= decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= decllow( max( eval( lhs ), eval( orhs ) ) );
            odres_  %= decllow( max( eval( lhs ), eval( orhs ) ) );
            sres_   %= decllow( max( eval( lhs ), eval( orhs ) ) );
            osres_  %= decllow( max( eval( lhs ), eval( orhs ) ) );
            refres_ %= decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= decllow( max( eval( olhs ), eval( rhs ) ) );
            odres_  %= decllow( max( eval( olhs ), eval( rhs ) ) );
            sres_   %= decllow( max( eval( olhs ), eval( rhs ) ) );
            osres_  %= decllow( max( eval( olhs ), eval( rhs ) ) );
            refres_ %= decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= decllow( max( eval( olhs ), eval( orhs ) ) );
            odres_  %= decllow( max( eval( olhs ), eval( orhs ) ) );
            sres_   %= decllow( max( eval( olhs ), eval( orhs ) ) );
            osres_  %= decllow( max( eval( olhs ), eval( orhs ) ) );
            refres_ %= decllow( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the lower dense matrix/dense matrix maximum operation.
//
// \return void
//
// This function is called in case the lower matrix/matrix maximum operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testDeclLowOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the upper dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the upper matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testDeclUppOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLUPP_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLUPP_OPERATION > 1 )
   {
      if( lhs_.rows() != lhs_.columns() )
         return;


      //=====================================================================================
      // Test-specific setup of the left-hand side operand
      //=====================================================================================

      MT1 lhs( lhs_ );

      blaze::resetLower( lhs );

      OMT1 olhs( lhs );


      //=====================================================================================
      // Test-specific setup of the right-hand side operand
      //=====================================================================================

      MT2 rhs( rhs_ );

      blaze::resetLower( rhs );

      OMT2 orhs( rhs );


      //=====================================================================================
      // Test-specific setup of the reference matrix
      //=====================================================================================

      RT ref( max( lhs, rhs ) );


      //=====================================================================================
      // Declupp maximum
      //=====================================================================================

      // Declupp maximum with the given matrices
      {
         test_  = "Declupp maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = declupp( max( lhs, rhs ) );
            odres_  = declupp( max( lhs, rhs ) );
            sres_   = declupp( max( lhs, rhs ) );
            osres_  = declupp( max( lhs, rhs ) );
            refres_ = declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declupp( max( lhs, orhs ) );
            odres_  = declupp( max( lhs, orhs ) );
            sres_   = declupp( max( lhs, orhs ) );
            osres_  = declupp( max( lhs, orhs ) );
            refres_ = declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declupp( max( olhs, rhs ) );
            odres_  = declupp( max( olhs, rhs ) );
            sres_   = declupp( max( olhs, rhs ) );
            osres_  = declupp( max( olhs, rhs ) );
            refres_ = declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declupp( max( olhs, orhs ) );
            odres_  = declupp( max( olhs, orhs ) );
            sres_   = declupp( max( olhs, orhs ) );
            osres_  = declupp( max( olhs, orhs ) );
            refres_ = declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declupp maximum with evaluated matrices
      {
         test_  = "Declupp maximum with evaluated left-hand side matrix";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = declupp( max( eval( lhs ), eval( rhs ) ) );
            odres_  = declupp( max( eval( lhs ), eval( rhs ) ) );
            sres_   = declupp( max( eval( lhs ), eval( rhs ) ) );
            osres_  = declupp( max( eval( lhs ), eval( rhs ) ) );
            refres_ = declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declupp( max( eval( lhs ), eval( orhs ) ) );
            odres_  = declupp( max( eval( lhs ), eval( orhs ) ) );
            sres_   = declupp( max( eval( lhs ), eval( orhs ) ) );
            osres_  = declupp( max( eval( lhs ), eval( orhs ) ) );
            refres_ = declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declupp( max( eval( olhs ), eval( rhs ) ) );
            odres_  = declupp( max( eval( olhs ), eval( rhs ) ) );
            sres_   = declupp( max( eval( olhs ), eval( rhs ) ) );
            osres_  = declupp( max( eval( olhs ), eval( rhs ) ) );
            refres_ = declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declupp( max( eval( olhs ), eval( orhs ) ) );
            odres_  = declupp( max( eval( olhs ), eval( orhs ) ) );
            sres_   = declupp( max( eval( olhs ), eval( orhs ) ) );
            osres_  = declupp( max( eval( olhs ), eval( orhs ) ) );
            refres_ = declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declupp maximum with addition assignment
      //=====================================================================================

      // Declupp maximum with addition assignment with the given matrices
      {
         test_  = "Declupp maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declupp( max( lhs, rhs ) );
            odres_  += declupp( max( lhs, rhs ) );
            sres_   += declupp( max( lhs, rhs ) );
            osres_  += declupp( max( lhs, rhs ) );
            refres_ += declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declupp( max( lhs, orhs ) );
            odres_  += declupp( max( lhs, orhs ) );
            sres_   += declupp( max( lhs, orhs ) );
            osres_  += declupp( max( lhs, orhs ) );
            refres_ += declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declupp( max( olhs, rhs ) );
            odres_  += declupp( max( olhs, rhs ) );
            sres_   += declupp( max( olhs, rhs ) );
            osres_  += declupp( max( olhs, rhs ) );
            refres_ += declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declupp( max( olhs, orhs ) );
            odres_  += declupp( max( olhs, orhs ) );
            sres_   += declupp( max( olhs, orhs ) );
            osres_  += declupp( max( olhs, orhs ) );
            refres_ += declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declupp maximum with addition assignment with evaluated matrices
      {
         test_  = "Declupp maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declupp( max( eval( lhs ), eval( rhs ) ) );
            odres_  += declupp( max( eval( lhs ), eval( rhs ) ) );
            sres_   += declupp( max( eval( lhs ), eval( rhs ) ) );
            osres_  += declupp( max( eval( lhs ), eval( rhs ) ) );
            refres_ += declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declupp( max( eval( lhs ), eval( orhs ) ) );
            odres_  += declupp( max( eval( lhs ), eval( orhs ) ) );
            sres_   += declupp( max( eval( lhs ), eval( orhs ) ) );
            osres_  += declupp( max( eval( lhs ), eval( orhs ) ) );
            refres_ += declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declupp( max( eval( olhs ), eval( rhs ) ) );
            odres_  += declupp( max( eval( olhs ), eval( rhs ) ) );
            sres_   += declupp( max( eval( olhs ), eval( rhs ) ) );
            osres_  += declupp( max( eval( olhs ), eval( rhs ) ) );
            refres_ += declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declupp( max( eval( olhs ), eval( orhs ) ) );
            odres_  += declupp( max( eval( olhs ), eval( orhs ) ) );
            sres_   += declupp( max( eval( olhs ), eval( orhs ) ) );
            osres_  += declupp( max( eval( olhs ), eval( orhs ) ) );
            refres_ += declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declupp maximum with subtraction assignment
      //=====================================================================================

      // Declupp maximum with subtraction assignment with the given matrices
      {
         test_  = "Declupp maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declupp( max( lhs, rhs ) );
            odres_  -= declupp( max( lhs, rhs ) );
            sres_   -= declupp( max( lhs, rhs ) );
            osres_  -= declupp( max( lhs, rhs ) );
            refres_ -= declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declupp( max( lhs, orhs ) );
            odres_  -= declupp( max( lhs, orhs ) );
            sres_   -= declupp( max( lhs, orhs ) );
            osres_  -= declupp( max( lhs, orhs ) );
            refres_ -= declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declupp( max( olhs, rhs ) );
            odres_  -= declupp( max( olhs, rhs ) );
            sres_   -= declupp( max( olhs, rhs ) );
            osres_  -= declupp( max( olhs, rhs ) );
            refres_ -= declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declupp( max( olhs, orhs ) );
            odres_  -= declupp( max( olhs, orhs ) );
            sres_   -= declupp( max( olhs, orhs ) );
            osres_  -= declupp( max( olhs, orhs ) );
            refres_ -= declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declupp maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Declupp maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declupp( max( eval( lhs ), eval( rhs ) ) );
            odres_  -= declupp( max( eval( lhs ), eval( rhs ) ) );
            sres_   -= declupp( max( eval( lhs ), eval( rhs ) ) );
            osres_  -= declupp( max( eval( lhs ), eval( rhs ) ) );
            refres_ -= declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declupp( max( eval( lhs ), eval( orhs ) ) );
            odres_  -= declupp( max( eval( lhs ), eval( orhs ) ) );
            sres_   -= declupp( max( eval( lhs ), eval( orhs ) ) );
            osres_  -= declupp( max( eval( lhs ), eval( orhs ) ) );
            refres_ -= declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declupp( max( eval( olhs ), eval( rhs ) ) );
            odres_  -= declupp( max( eval( olhs ), eval( rhs ) ) );
            sres_   -= declupp( max( eval( olhs ), eval( rhs ) ) );
            osres_  -= declupp( max( eval( olhs ), eval( rhs ) ) );
            refres_ -= declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declupp( max( eval( olhs ), eval( orhs ) ) );
            odres_  -= declupp( max( eval( olhs ), eval( orhs ) ) );
            sres_   -= declupp( max( eval( olhs ), eval( orhs ) ) );
            osres_  -= declupp( max( eval( olhs ), eval( orhs ) ) );
            refres_ -= declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declupp maximum with Schur product assignment
      //=====================================================================================

      // Declupp maximum with Schur product assignment with the given matrices
      {
         test_  = "Declupp maximum with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declupp( max( lhs, rhs ) );
            odres_  %= declupp( max( lhs, rhs ) );
            sres_   %= declupp( max( lhs, rhs ) );
            osres_  %= declupp( max( lhs, rhs ) );
            refres_ %= declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declupp( max( lhs, orhs ) );
            odres_  %= declupp( max( lhs, orhs ) );
            sres_   %= declupp( max( lhs, orhs ) );
            osres_  %= declupp( max( lhs, orhs ) );
            refres_ %= declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declupp( max( olhs, rhs ) );
            odres_  %= declupp( max( olhs, rhs ) );
            sres_   %= declupp( max( olhs, rhs ) );
            osres_  %= declupp( max( olhs, rhs ) );
            refres_ %= declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declupp( max( olhs, orhs ) );
            odres_  %= declupp( max( olhs, orhs ) );
            sres_   %= declupp( max( olhs, orhs ) );
            osres_  %= declupp( max( olhs, orhs ) );
            refres_ %= declupp( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declupp maximum with Schur product assignment with evaluated matrices
      {
         test_  = "Declupp maximum with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declupp( max( eval( lhs ), eval( rhs ) ) );
            odres_  %= declupp( max( eval( lhs ), eval( rhs ) ) );
            sres_   %= declupp( max( eval( lhs ), eval( rhs ) ) );
            osres_  %= declupp( max( eval( lhs ), eval( rhs ) ) );
            refres_ %= declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declupp( max( eval( lhs ), eval( orhs ) ) );
            odres_  %= declupp( max( eval( lhs ), eval( orhs ) ) );
            sres_   %= declupp( max( eval( lhs ), eval( orhs ) ) );
            osres_  %= declupp( max( eval( lhs ), eval( orhs ) ) );
            refres_ %= declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declupp( max( eval( olhs ), eval( rhs ) ) );
            odres_  %= declupp( max( eval( olhs ), eval( rhs ) ) );
            sres_   %= declupp( max( eval( olhs ), eval( rhs ) ) );
            osres_  %= declupp( max( eval( olhs ), eval( rhs ) ) );
            refres_ %= declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declupp( max( eval( olhs ), eval( orhs ) ) );
            odres_  %= declupp( max( eval( olhs ), eval( orhs ) ) );
            sres_   %= declupp( max( eval( olhs ), eval( orhs ) ) );
            osres_  %= declupp( max( eval( olhs ), eval( orhs ) ) );
            refres_ %= declupp( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the upper dense matrix/dense matrix maximum operation.
//
// \return void
//
// This function is called in case the upper matrix/matrix maximum operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testDeclUppOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the diagonal dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the diagonal matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testDeclDiagOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLDIAG_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLDIAG_OPERATION > 1 )
   {
      if( lhs_.rows() != lhs_.columns() )
         return;


      //=====================================================================================
      // Test-specific setup of the left-hand side operand
      //=====================================================================================

      MT1 lhs( lhs_ );

      blaze::resetLower( lhs );
      blaze::resetUpper( lhs );

      OMT1 olhs( lhs );


      //=====================================================================================
      // Test-specific setup of the right-hand side operand
      //=====================================================================================

      MT2 rhs( rhs_ );

      blaze::resetLower( rhs );
      blaze::resetUpper( rhs );

      OMT2 orhs( rhs );


      //=====================================================================================
      // Test-specific setup of the reference matrix
      //=====================================================================================

      RT ref( max( lhs, rhs ) );


      //=====================================================================================
      // Decldiag maximum
      //=====================================================================================

      // Decldiag maximum with the given matrices
      {
         test_  = "Decldiag maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = decldiag( max( lhs, rhs ) );
            odres_  = decldiag( max( lhs, rhs ) );
            sres_   = decldiag( max( lhs, rhs ) );
            osres_  = decldiag( max( lhs, rhs ) );
            refres_ = decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = decldiag( max( lhs, orhs ) );
            odres_  = decldiag( max( lhs, orhs ) );
            sres_   = decldiag( max( lhs, orhs ) );
            osres_  = decldiag( max( lhs, orhs ) );
            refres_ = decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = decldiag( max( olhs, rhs ) );
            odres_  = decldiag( max( olhs, rhs ) );
            sres_   = decldiag( max( olhs, rhs ) );
            osres_  = decldiag( max( olhs, rhs ) );
            refres_ = decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = decldiag( max( olhs, orhs ) );
            odres_  = decldiag( max( olhs, orhs ) );
            sres_   = decldiag( max( olhs, orhs ) );
            osres_  = decldiag( max( olhs, orhs ) );
            refres_ = decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decldiag maximum with evaluated matrices
      {
         test_  = "Decldiag maximum with evaluated left-hand side matrix";
         error_ = "Failed maximum operation";

         try {
            initResults();
            dres_   = decldiag( max( eval( lhs ), eval( rhs ) ) );
            odres_  = decldiag( max( eval( lhs ), eval( rhs ) ) );
            sres_   = decldiag( max( eval( lhs ), eval( rhs ) ) );
            osres_  = decldiag( max( eval( lhs ), eval( rhs ) ) );
            refres_ = decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = decldiag( max( eval( lhs ), eval( orhs ) ) );
            odres_  = decldiag( max( eval( lhs ), eval( orhs ) ) );
            sres_   = decldiag( max( eval( lhs ), eval( orhs ) ) );
            osres_  = decldiag( max( eval( lhs ), eval( orhs ) ) );
            refres_ = decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = decldiag( max( eval( olhs ), eval( rhs ) ) );
            odres_  = decldiag( max( eval( olhs ), eval( rhs ) ) );
            sres_   = decldiag( max( eval( olhs ), eval( rhs ) ) );
            osres_  = decldiag( max( eval( olhs ), eval( rhs ) ) );
            refres_ = decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = decldiag( max( eval( olhs ), eval( orhs ) ) );
            odres_  = decldiag( max( eval( olhs ), eval( orhs ) ) );
            sres_   = decldiag( max( eval( olhs ), eval( orhs ) ) );
            osres_  = decldiag( max( eval( olhs ), eval( orhs ) ) );
            refres_ = decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Decldiag maximum with addition assignment
      //=====================================================================================

      // Decldiag maximum with addition assignment with the given matrices
      {
         test_  = "Decldiag maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += decldiag( max( lhs, rhs ) );
            odres_  += decldiag( max( lhs, rhs ) );
            sres_   += decldiag( max( lhs, rhs ) );
            osres_  += decldiag( max( lhs, rhs ) );
            refres_ += decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += decldiag( max( lhs, orhs ) );
            odres_  += decldiag( max( lhs, orhs ) );
            sres_   += decldiag( max( lhs, orhs ) );
            osres_  += decldiag( max( lhs, orhs ) );
            refres_ += decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += decldiag( max( olhs, rhs ) );
            odres_  += decldiag( max( olhs, rhs ) );
            sres_   += decldiag( max( olhs, rhs ) );
            osres_  += decldiag( max( olhs, rhs ) );
            refres_ += decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += decldiag( max( olhs, orhs ) );
            odres_  += decldiag( max( olhs, orhs ) );
            sres_   += decldiag( max( olhs, orhs ) );
            osres_  += decldiag( max( olhs, orhs ) );
            refres_ += decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decldiag maximum with addition assignment with evaluated matrices
      {
         test_  = "Decldiag maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += decldiag( max( eval( lhs ), eval( rhs ) ) );
            odres_  += decldiag( max( eval( lhs ), eval( rhs ) ) );
            sres_   += decldiag( max( eval( lhs ), eval( rhs ) ) );
            osres_  += decldiag( max( eval( lhs ), eval( rhs ) ) );
            refres_ += decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += decldiag( max( eval( lhs ), eval( orhs ) ) );
            odres_  += decldiag( max( eval( lhs ), eval( orhs ) ) );
            sres_   += decldiag( max( eval( lhs ), eval( orhs ) ) );
            osres_  += decldiag( max( eval( lhs ), eval( orhs ) ) );
            refres_ += decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += decldiag( max( eval( olhs ), eval( rhs ) ) );
            odres_  += decldiag( max( eval( olhs ), eval( rhs ) ) );
            sres_   += decldiag( max( eval( olhs ), eval( rhs ) ) );
            osres_  += decldiag( max( eval( olhs ), eval( rhs ) ) );
            refres_ += decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += decldiag( max( eval( olhs ), eval( orhs ) ) );
            odres_  += decldiag( max( eval( olhs ), eval( orhs ) ) );
            sres_   += decldiag( max( eval( olhs ), eval( orhs ) ) );
            osres_  += decldiag( max( eval( olhs ), eval( orhs ) ) );
            refres_ += decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Decldiag maximum with subtraction assignment
      //=====================================================================================

      // Decldiag maximum with subtraction assignment with the given matrices
      {
         test_  = "Decldiag maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= decldiag( max( lhs, rhs ) );
            odres_  -= decldiag( max( lhs, rhs ) );
            sres_   -= decldiag( max( lhs, rhs ) );
            osres_  -= decldiag( max( lhs, rhs ) );
            refres_ -= decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= decldiag( max( lhs, orhs ) );
            odres_  -= decldiag( max( lhs, orhs ) );
            sres_   -= decldiag( max( lhs, orhs ) );
            osres_  -= decldiag( max( lhs, orhs ) );
            refres_ -= decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= decldiag( max( olhs, rhs ) );
            odres_  -= decldiag( max( olhs, rhs ) );
            sres_   -= decldiag( max( olhs, rhs ) );
            osres_  -= decldiag( max( olhs, rhs ) );
            refres_ -= decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= decldiag( max( olhs, orhs ) );
            odres_  -= decldiag( max( olhs, orhs ) );
            sres_   -= decldiag( max( olhs, orhs ) );
            osres_  -= decldiag( max( olhs, orhs ) );
            refres_ -= decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decldiag maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Decldiag maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= decldiag( max( eval( lhs ), eval( rhs ) ) );
            odres_  -= decldiag( max( eval( lhs ), eval( rhs ) ) );
            sres_   -= decldiag( max( eval( lhs ), eval( rhs ) ) );
            osres_  -= decldiag( max( eval( lhs ), eval( rhs ) ) );
            refres_ -= decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= decldiag( max( eval( lhs ), eval( orhs ) ) );
            odres_  -= decldiag( max( eval( lhs ), eval( orhs ) ) );
            sres_   -= decldiag( max( eval( lhs ), eval( orhs ) ) );
            osres_  -= decldiag( max( eval( lhs ), eval( orhs ) ) );
            refres_ -= decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= decldiag( max( eval( olhs ), eval( rhs ) ) );
            odres_  -= decldiag( max( eval( olhs ), eval( rhs ) ) );
            sres_   -= decldiag( max( eval( olhs ), eval( rhs ) ) );
            osres_  -= decldiag( max( eval( olhs ), eval( rhs ) ) );
            refres_ -= decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= decldiag( max( eval( olhs ), eval( orhs ) ) );
            odres_  -= decldiag( max( eval( olhs ), eval( orhs ) ) );
            sres_   -= decldiag( max( eval( olhs ), eval( orhs ) ) );
            osres_  -= decldiag( max( eval( olhs ), eval( orhs ) ) );
            refres_ -= decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Decldiag maximum with Schur product assignment
      //=====================================================================================

      // Decldiag maximum with Schur product assignment with the given matrices
      {
         test_  = "Decldiag maximum with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= decldiag( max( lhs, rhs ) );
            odres_  %= decldiag( max( lhs, rhs ) );
            sres_   %= decldiag( max( lhs, rhs ) );
            osres_  %= decldiag( max( lhs, rhs ) );
            refres_ %= decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= decldiag( max( lhs, orhs ) );
            odres_  %= decldiag( max( lhs, orhs ) );
            sres_   %= decldiag( max( lhs, orhs ) );
            osres_  %= decldiag( max( lhs, orhs ) );
            refres_ %= decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= decldiag( max( olhs, rhs ) );
            odres_  %= decldiag( max( olhs, rhs ) );
            sres_   %= decldiag( max( olhs, rhs ) );
            osres_  %= decldiag( max( olhs, rhs ) );
            refres_ %= decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= decldiag( max( olhs, orhs ) );
            odres_  %= decldiag( max( olhs, orhs ) );
            sres_   %= decldiag( max( olhs, orhs ) );
            osres_  %= decldiag( max( olhs, orhs ) );
            refres_ %= decldiag( ref );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decldiag maximum with Schur product assignment with evaluated matrices
      {
         test_  = "Decldiag maximum with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= decldiag( max( eval( lhs ), eval( rhs ) ) );
            odres_  %= decldiag( max( eval( lhs ), eval( rhs ) ) );
            sres_   %= decldiag( max( eval( lhs ), eval( rhs ) ) );
            osres_  %= decldiag( max( eval( lhs ), eval( rhs ) ) );
            refres_ %= decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= decldiag( max( eval( lhs ), eval( orhs ) ) );
            odres_  %= decldiag( max( eval( lhs ), eval( orhs ) ) );
            sres_   %= decldiag( max( eval( lhs ), eval( orhs ) ) );
            osres_  %= decldiag( max( eval( lhs ), eval( orhs ) ) );
            refres_ %= decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= decldiag( max( eval( olhs ), eval( rhs ) ) );
            odres_  %= decldiag( max( eval( olhs ), eval( rhs ) ) );
            sres_   %= decldiag( max( eval( olhs ), eval( rhs ) ) );
            osres_  %= decldiag( max( eval( olhs ), eval( rhs ) ) );
            refres_ %= decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= decldiag( max( eval( olhs ), eval( orhs ) ) );
            odres_  %= decldiag( max( eval( olhs ), eval( orhs ) ) );
            sres_   %= decldiag( max( eval( olhs ), eval( orhs ) ) );
            osres_  %= decldiag( max( eval( olhs ), eval( orhs ) ) );
            refres_ %= decldiag( eval( ref ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the diagonal dense matrix/dense matrix maximum operation.
//
// \return void
//
// This function is called in case the diagonal matrix/matrix maximum operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testDeclDiagOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the submatrix-wise dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the submatrix-wise matrix maximum with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the maximum operation or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testSubmatrixOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION > 1 )
   {
      if( lhs_.rows() == 0UL || lhs_.columns() == 0UL )
         return;


      //=====================================================================================
      // Submatrix-wise maximum
      //=====================================================================================

      // Submatrix-wise maximum with the given matrices
      {
         test_  = "Submatrix-wise maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( ref_             , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( ref_              , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( ref_              , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( ref_               , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Submatrix-wise maximum with evaluated matrices
      {
         test_  = "Submatrix-wise maximum with evaluated matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( ref_ )                     , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( ref_ )                      , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( ref_ )                      , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( ref_ )                       , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Submatrix-wise maximum with addition assignment
      //=====================================================================================

      // Submatrix-wise maximum with addition assignment with the given matrices
      {
         test_  = "Submatrix-wise maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( ref_             , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( ref_              , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( ref_              , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( ref_               , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Submatrix-wise maximum with addition assignment with evaluated matrices
      {
         test_  = "Submatrix-wise maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( ref_ )                     , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( ref_ )                      , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( ref_ )                      , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( ref_ )                       , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Submatrix-wise maximum with subtraction assignment
      //=====================================================================================

      // Submatrix-wise maximum with subtraction assignment with the given matrices
      {
         test_  = "Submatrix-wise maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( ref_             , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( ref_              , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( ref_              , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( ref_               , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Submatrix-wise maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Submatrix-wise maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( ref_ )                     , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( ref_ )                      , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( ref_ )                      , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( ref_ )                       , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Submatrix-wise maximum with Schur product assignment
      //=====================================================================================

      // Submatrix-wise maximum with Schur product assignment with the given matrices
      {
         test_  = "Submatrix-wise maximum with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( max( lhs_, rhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( ref_             , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( max( lhs_, orhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( ref_              , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( max( olhs_, rhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( ref_              , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( max( olhs_, orhs_ ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( ref_               , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Submatrix-wise maximum with Schur product assignment with evaluated matrices
      {
         test_  = "Submatrix-wise maximum with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( max( eval( lhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( eval( ref_ )                     , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( max( eval( lhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( eval( ref_ )                      , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( max( eval( olhs_ ), eval( rhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( eval( ref_ )                      , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<orhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, orhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( max( eval( olhs_ ), eval( orhs_ ) ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( eval( ref_ )                       , row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds submatrix construction (invalid number of rows)
      {
         test_  = "Out-of-bounds submatrix construction (invalid number of rows)";
         error_ = "Setup of out-of-bounds submatrix succeeded";

         try {
            auto sm = submatrix( max( lhs_, rhs_ ), 1UL, 0UL, lhs_.rows(), lhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( max( lhs_, orhs_ ), 1UL, 0UL, lhs_.rows(), lhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( max( olhs_, rhs_ ), 1UL, 0UL, olhs_.rows(), olhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( max( olhs_, orhs_ ), 1UL, 0UL, olhs_.rows(), olhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
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
            auto sm = submatrix( max( lhs_, rhs_ ), 0UL, 1UL, lhs_.rows(), lhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( max( lhs_, orhs_ ), 0UL, 1UL, lhs_.rows(), lhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( max( olhs_, rhs_ ), 0UL, 1UL, olhs_.rows(), olhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( max( olhs_, orhs_ ), 0UL, 1UL, olhs_.rows(), olhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
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
            auto sm = submatrix( max( lhs_, rhs_ ), lhs_.rows(), 0UL, 1UL, lhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( max( lhs_, orhs_ ), lhs_.rows(), 0UL, 1UL, lhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( max( olhs_, rhs_ ), olhs_.rows(), 0UL, 1UL, olhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( max( olhs_, orhs_ ), olhs_.rows(), 0UL, 1UL, olhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
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
            auto sm = submatrix( max( lhs_, rhs_ ), 0UL, lhs_.columns(), lhs_.rows(), 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( max( lhs_, orhs_ ), 0UL, lhs_.columns(), lhs_.rows(), 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( max( olhs_, rhs_ ), 0UL, olhs_.columns(), olhs_.rows(), 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( max( olhs_, orhs_ ), 0UL, olhs_.columns(), olhs_.rows(), 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
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
/*!\brief Skipping the submatrix-wise dense matrix/dense matrix maximum operation.
//
// \return void
//
// This function is called in case the submatrix-wise matrix/matrix maximum operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testSubmatrixOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the row-wise dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the row-wise matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testRowOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_ROW_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ROW_OPERATION > 1 )
   {
      if( lhs_.rows() == 0UL )
         return;


      //=====================================================================================
      // Row-wise maximum
      //=====================================================================================

      // Row-wise maximum with the given matrices
      {
         test_  = "Row-wise maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( max( lhs_, rhs_ ), i );
               row( odres_ , i ) = row( max( lhs_, rhs_ ), i );
               row( sres_  , i ) = row( max( lhs_, rhs_ ), i );
               row( osres_ , i ) = row( max( lhs_, rhs_ ), i );
               row( refres_, i ) = row( ref_             , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( max( lhs_, orhs_ ), i );
               row( odres_ , i ) = row( max( lhs_, orhs_ ), i );
               row( sres_  , i ) = row( max( lhs_, orhs_ ), i );
               row( osres_ , i ) = row( max( lhs_, orhs_ ), i );
               row( refres_, i ) = row( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( max( olhs_, rhs_ ), i );
               row( odres_ , i ) = row( max( olhs_, rhs_ ), i );
               row( sres_  , i ) = row( max( olhs_, rhs_ ), i );
               row( osres_ , i ) = row( max( olhs_, rhs_ ), i );
               row( refres_, i ) = row( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( max( olhs_, orhs_ ), i );
               row( odres_ , i ) = row( max( olhs_, orhs_ ), i );
               row( sres_  , i ) = row( max( olhs_, orhs_ ), i );
               row( osres_ , i ) = row( max( olhs_, orhs_ ), i );
               row( refres_, i ) = row( ref_               , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Row-wise maximum with evaluated matrices
      {
         test_  = "Row-wise maximum with evaluated matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( odres_ , i ) = row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( sres_  , i ) = row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( osres_ , i ) = row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( refres_, i ) = row( eval( ref_ )                     , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( odres_ , i ) = row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( sres_  , i ) = row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( osres_ , i ) = row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( refres_, i ) = row( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( odres_ , i ) = row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( sres_  , i ) = row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( osres_ , i ) = row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( refres_, i ) = row( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( odres_ , i ) = row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( sres_  , i ) = row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( osres_ , i ) = row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( refres_, i ) = row( eval( ref_ )                       , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Row-wise maximum with addition assignment
      //=====================================================================================

      // Row-wise maximum with addition assignment with the given matrices
      {
         test_  = "Row-wise maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( max( lhs_, rhs_ ), i );
               row( odres_ , i ) += row( max( lhs_, rhs_ ), i );
               row( sres_  , i ) += row( max( lhs_, rhs_ ), i );
               row( osres_ , i ) += row( max( lhs_, rhs_ ), i );
               row( refres_, i ) += row( ref_             , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( max( lhs_, orhs_ ), i );
               row( odres_ , i ) += row( max( lhs_, orhs_ ), i );
               row( sres_  , i ) += row( max( lhs_, orhs_ ), i );
               row( osres_ , i ) += row( max( lhs_, orhs_ ), i );
               row( refres_, i ) += row( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( max( olhs_, rhs_ ), i );
               row( odres_ , i ) += row( max( olhs_, rhs_ ), i );
               row( sres_  , i ) += row( max( olhs_, rhs_ ), i );
               row( osres_ , i ) += row( max( olhs_, rhs_ ), i );
               row( refres_, i ) += row( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( max( olhs_, orhs_ ), i );
               row( odres_ , i ) += row( max( olhs_, orhs_ ), i );
               row( sres_  , i ) += row( max( olhs_, orhs_ ), i );
               row( osres_ , i ) += row( max( olhs_, orhs_ ), i );
               row( refres_, i ) += row( ref_               , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Row-wise maximum with addition assignment with evaluated matrices
      {
         test_  = "Row-wise maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( odres_ , i ) += row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( sres_  , i ) += row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( osres_ , i ) += row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( refres_, i ) += row( eval( ref_ )                     , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( odres_ , i ) += row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( sres_  , i ) += row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( osres_ , i ) += row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( refres_, i ) += row( eval( ref_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( odres_ , i ) += row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( sres_  , i ) += row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( osres_ , i ) += row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( refres_, i ) += row( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( odres_ , i ) += row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( sres_  , i ) += row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( osres_ , i ) += row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( refres_, i ) += row( eval( ref_ )                       , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Row-wise maximum with subtraction assignment
      //=====================================================================================

      // Row-wise maximum with subtraction assignment with the given matrices
      {
         test_  = "Row-wise maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( max( lhs_, rhs_ ), i );
               row( odres_ , i ) -= row( max( lhs_, rhs_ ), i );
               row( sres_  , i ) -= row( max( lhs_, rhs_ ), i );
               row( osres_ , i ) -= row( max( lhs_, rhs_ ), i );
               row( refres_, i ) -= row( ref_             , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( max( lhs_, orhs_ ), i );
               row( odres_ , i ) -= row( max( lhs_, orhs_ ), i );
               row( sres_  , i ) -= row( max( lhs_, orhs_ ), i );
               row( osres_ , i ) -= row( max( lhs_, orhs_ ), i );
               row( refres_, i ) -= row( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( max( olhs_, rhs_ ), i );
               row( odres_ , i ) -= row( max( olhs_, rhs_ ), i );
               row( sres_  , i ) -= row( max( olhs_, rhs_ ), i );
               row( osres_ , i ) -= row( max( olhs_, rhs_ ), i );
               row( refres_, i ) -= row( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( max( olhs_, orhs_ ), i );
               row( odres_ , i ) -= row( max( olhs_, orhs_ ), i );
               row( sres_  , i ) -= row( max( olhs_, orhs_ ), i );
               row( osres_ , i ) -= row( max( olhs_, orhs_ ), i );
               row( refres_, i ) -= row( ref_               , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Row-wise maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Row-wise maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( odres_ , i ) -= row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( sres_  , i ) -= row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( osres_ , i ) -= row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( refres_, i ) -= row( eval( ref_ )                     , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( odres_ , i ) -= row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( sres_  , i ) -= row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( osres_ , i ) -= row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( refres_, i ) -= row( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( odres_ , i ) -= row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( sres_  , i ) -= row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( osres_ , i ) -= row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( refres_, i ) -= row( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( odres_ , i ) -= row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( sres_  , i ) -= row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( osres_ , i ) -= row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( refres_, i ) -= row( eval( ref_ )                       , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Row-wise maximum with multiplication assignment
      //=====================================================================================

      // Row-wise maximum with multiplication assignment with the given matrices
      {
         test_  = "Row-wise maximum with multiplication assignment with the given matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( max( lhs_, rhs_ ), i );
               row( odres_ , i ) *= row( max( lhs_, rhs_ ), i );
               row( sres_  , i ) *= row( max( lhs_, rhs_ ), i );
               row( osres_ , i ) *= row( max( lhs_, rhs_ ), i );
               row( refres_, i ) *= row( ref_             , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( max( lhs_, orhs_ ), i );
               row( odres_ , i ) *= row( max( lhs_, orhs_ ), i );
               row( sres_  , i ) *= row( max( lhs_, orhs_ ), i );
               row( osres_ , i ) *= row( max( lhs_, orhs_ ), i );
               row( refres_, i ) *= row( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( max( olhs_, rhs_ ), i );
               row( odres_ , i ) *= row( max( olhs_, rhs_ ), i );
               row( sres_  , i ) *= row( max( olhs_, rhs_ ), i );
               row( osres_ , i ) *= row( max( olhs_, rhs_ ), i );
               row( refres_, i ) *= row( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( max( olhs_, orhs_ ), i );
               row( odres_ , i ) *= row( max( olhs_, orhs_ ), i );
               row( sres_  , i ) *= row( max( olhs_, orhs_ ), i );
               row( osres_ , i ) *= row( max( olhs_, orhs_ ), i );
               row( refres_, i ) *= row( ref_               , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Row-wise maximum with multiplication assignment with evaluated matrices
      {
         test_  = "Row-wise maximum with multiplication assignment with evaluated matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( odres_ , i ) *= row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( sres_  , i ) *= row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( osres_ , i ) *= row( max( eval( lhs_ ), eval( rhs_ ) ), i );
               row( refres_, i ) *= row( eval( ref_ )                     , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( odres_ , i ) *= row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( sres_  , i ) *= row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( osres_ , i ) *= row( max( eval( lhs_ ), eval( orhs_ ) ), i );
               row( refres_, i ) *= row( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( odres_ , i ) *= row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( sres_  , i ) *= row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( osres_ , i ) *= row( max( eval( olhs_ ), eval( rhs_ ) ), i );
               row( refres_, i ) *= row( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( odres_ , i ) *= row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( sres_  , i ) *= row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( osres_ , i ) *= row( max( eval( olhs_ ), eval( orhs_ ) ), i );
               row( refres_, i ) *= row( eval( ref_ )                       , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid row index)
      {
         test_  = "Out-of-bounds row construction (invalid row index)";
         error_ = "Setup of out-of-bounds row succeeded";

         try {
            auto r = row( max( lhs_, rhs_ ), lhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = row( max( lhs_, orhs_ ), lhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = row( max( olhs_, rhs_ ), olhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = row( max( olhs_, orhs_ ), olhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
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
/*!\brief Skipping the row-wise dense matrix/dense matrix maximum operation.
//
// \return void
//
// This function is called in case the row-wise matrix/matrix maximum operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testRowOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the rows-wise dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the rows-wise matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testRowsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_ROWS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ROWS_OPERATION > 1 )
   {
      if( lhs_.rows() == 0UL )
         return;


      std::vector<size_t> indices( lhs_.rows() );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );


      //=====================================================================================
      // Rows-wise maximum
      //=====================================================================================

      // Rows-wise maximum with the given matrices
      {
         test_  = "Rows-wise maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Rows-wise maximum with evaluated matrices
      {
         test_  = "Rows-wise maximum with evaluated matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Rows-wise maximum with addition assignment
      //=====================================================================================

      // Rows-wise maximum with addition assignment with the given matrices
      {
         test_  = "Rows-wise maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Rows-wise maximum with addition assignment with evaluated matrices
      {
         test_  = "Rows-wise maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Rows-wise maximum with subtraction assignment
      //=====================================================================================

      // Rows-wise maximum with subtraction assignment with the given matrices
      {
         test_  = "Rows-wise maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Rows-wise maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Rows-wise maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Rows-wise maximum with Schur product assignment
      //=====================================================================================

      // Rows-wise maximum with Schur product assignment with the given matrices
      {
         test_  = "Rows-wise maximum with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( max( lhs_, rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( max( lhs_, orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( max( olhs_, rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( max( olhs_, orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Rows-wise maximum with Schur product assignment with evaluated matrices
      {
         test_  = "Rows-wise maximum with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid row index; initializer_list)
      {
         test_  = "Out-of-bounds row selection construction (invalid row index; initializer_list)";
         error_ = "Setup of out-of-bounds row selection succeeded";

         try {
            auto r = rows( max( lhs_, rhs_ ), { lhs_.rows() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( max( lhs_, orhs_ ), { lhs_.rows() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( max( olhs_, rhs_ ), { olhs_.rows() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( max( olhs_, orhs_ ), { olhs_.rows() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
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
            auto r = rows( max( lhs_, rhs_ ), [index=lhs_.rows()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( max( lhs_, orhs_ ), [index=lhs_.rows()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( max( olhs_, rhs_ ), [index=olhs_.rows()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( max( olhs_, orhs_ ), [index=olhs_.rows()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
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
/*!\brief Skipping the rows-wise dense matrix/dense matrix maximum operation.
//
// \return void
//
// This function is called in case the rows-wise matrix/matrix maximum operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testRowsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the column-wise dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the column-wise matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testColumnOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION
   if( BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION > 1 )
   {
      if( lhs_.columns() == 0UL )
         return;


      //=====================================================================================
      // Column-wise maximum
      //=====================================================================================

      // Column-wise maximum with the given matrices
      {
         test_  = "Column-wise maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( max( lhs_, rhs_ ), j );
               column( odres_ , j ) = column( max( lhs_, rhs_ ), j );
               column( sres_  , j ) = column( max( lhs_, rhs_ ), j );
               column( osres_ , j ) = column( max( lhs_, rhs_ ), j );
               column( refres_, j ) = column( ref_             , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( max( lhs_, orhs_ ), j );
               column( odres_ , j ) = column( max( lhs_, orhs_ ), j );
               column( sres_  , j ) = column( max( lhs_, orhs_ ), j );
               column( osres_ , j ) = column( max( lhs_, orhs_ ), j );
               column( refres_, j ) = column( ref_              , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( max( olhs_, rhs_ ), j );
               column( odres_ , j ) = column( max( olhs_, rhs_ ), j );
               column( sres_  , j ) = column( max( olhs_, rhs_ ), j );
               column( osres_ , j ) = column( max( olhs_, rhs_ ), j );
               column( refres_, j ) = column( ref_              , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( max( olhs_, orhs_ ), j );
               column( odres_ , j ) = column( max( olhs_, orhs_ ), j );
               column( sres_  , j ) = column( max( olhs_, orhs_ ), j );
               column( osres_ , j ) = column( max( olhs_, orhs_ ), j );
               column( refres_, j ) = column( ref_               , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Column-wise maximum with evaluated matrices
      {
         test_  = "Column-wise maximum with evaluated matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( odres_ , j ) = column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( sres_  , j ) = column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( osres_ , j ) = column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( refres_, j ) = column( eval( ref_ )                     , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( odres_ , j ) = column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( sres_  , j ) = column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( osres_ , j ) = column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( refres_, j ) = column( eval( ref_ )                      , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( odres_ , j ) = column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( sres_  , j ) = column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( osres_ , j ) = column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( refres_, j ) = column( eval( ref_ )                      , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) = column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( odres_ , j ) = column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( sres_  , j ) = column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( osres_ , j ) = column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( refres_, j ) = column( eval( ref_ )                       , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Column-wise maximum with addition assignment
      //=====================================================================================

      // Column-wise maximum with addition assignment with the given matrices
      {
         test_  = "Column-wise maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( max( lhs_, rhs_ ), j );
               column( odres_ , j ) += column( max( lhs_, rhs_ ), j );
               column( sres_  , j ) += column( max( lhs_, rhs_ ), j );
               column( osres_ , j ) += column( max( lhs_, rhs_ ), j );
               column( refres_, j ) += column( ref_             , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( max( lhs_, orhs_ ), j );
               column( odres_ , j ) += column( max( lhs_, orhs_ ), j );
               column( sres_  , j ) += column( max( lhs_, orhs_ ), j );
               column( osres_ , j ) += column( max( lhs_, orhs_ ), j );
               column( refres_, j ) += column( ref_              , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( max( olhs_, rhs_ ), j );
               column( odres_ , j ) += column( max( olhs_, rhs_ ), j );
               column( sres_  , j ) += column( max( olhs_, rhs_ ), j );
               column( osres_ , j ) += column( max( olhs_, rhs_ ), j );
               column( refres_, j ) += column( ref_              , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( max( olhs_, orhs_ ), j );
               column( odres_ , j ) += column( max( olhs_, orhs_ ), j );
               column( sres_  , j ) += column( max( olhs_, orhs_ ), j );
               column( osres_ , j ) += column( max( olhs_, orhs_ ), j );
               column( refres_, j ) += column( ref_               , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Column-wise maximum with addition assignment with evaluated matrices
      {
         test_  = "Column-wise maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( odres_ , j ) += column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( sres_  , j ) += column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( osres_ , j ) += column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( refres_, j ) += column( eval( ref_ )                     , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( odres_ , j ) += column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( sres_  , j ) += column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( osres_ , j ) += column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( refres_, j ) += column( eval( ref_ )                      , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( odres_ , j ) += column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( sres_  , j ) += column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( osres_ , j ) += column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( refres_, j ) += column( eval( ref_ )                      , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) += column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( odres_ , j ) += column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( sres_  , j ) += column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( osres_ , j ) += column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( refres_, j ) += column( eval( ref_ )                       , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Column-wise maximum with subtraction assignment
      //=====================================================================================

      // Column-wise maximum with subtraction assignment with the given matrices
      {
         test_  = "Column-wise maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( max( lhs_, rhs_ ), j );
               column( odres_ , j ) -= column( max( lhs_, rhs_ ), j );
               column( sres_  , j ) -= column( max( lhs_, rhs_ ), j );
               column( osres_ , j ) -= column( max( lhs_, rhs_ ), j );
               column( refres_, j ) -= column( ref_             , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( max( lhs_, orhs_ ), j );
               column( odres_ , j ) -= column( max( lhs_, orhs_ ), j );
               column( sres_  , j ) -= column( max( lhs_, orhs_ ), j );
               column( osres_ , j ) -= column( max( lhs_, orhs_ ), j );
               column( refres_, j ) -= column( ref_              , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( max( olhs_, rhs_ ), j );
               column( odres_ , j ) -= column( max( olhs_, rhs_ ), j );
               column( sres_  , j ) -= column( max( olhs_, rhs_ ), j );
               column( osres_ , j ) -= column( max( olhs_, rhs_ ), j );
               column( refres_, j ) -= column( ref_              , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( max( olhs_, orhs_ ), j );
               column( odres_ , j ) -= column( max( olhs_, orhs_ ), j );
               column( sres_  , j ) -= column( max( olhs_, orhs_ ), j );
               column( osres_ , j ) -= column( max( olhs_, orhs_ ), j );
               column( refres_, j ) -= column( ref_               , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Column-wise maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Column-wise maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( odres_ , j ) -= column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( sres_  , j ) -= column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( osres_ , j ) -= column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( refres_, j ) -= column( eval( ref_ )                     , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( odres_ , j ) -= column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( sres_  , j ) -= column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( osres_ , j ) -= column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( refres_, j ) -= column( eval( ref_ )                      , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( odres_ , j ) -= column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( sres_  , j ) -= column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( osres_ , j ) -= column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( refres_, j ) -= column( eval( ref_ )                      , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( odres_ , j ) -= column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( sres_  , j ) -= column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( osres_ , j ) -= column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( refres_, j ) -= column( eval( ref_ )                       , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Column-wise maximum with multiplication assignment
      //=====================================================================================

      // Column-wise maximum with multiplication assignment with the given matrices
      {
         test_  = "Column-wise maximum with multiplication assignment with the given matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( max( lhs_, rhs_ ), j );
               column( odres_ , j ) *= column( max( lhs_, rhs_ ), j );
               column( sres_  , j ) *= column( max( lhs_, rhs_ ), j );
               column( osres_ , j ) *= column( max( lhs_, rhs_ ), j );
               column( refres_, j ) *= column( ref_             , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( max( lhs_, orhs_ ), j );
               column( odres_ , j ) *= column( max( lhs_, orhs_ ), j );
               column( sres_  , j ) *= column( max( lhs_, orhs_ ), j );
               column( osres_ , j ) *= column( max( lhs_, orhs_ ), j );
               column( refres_, j ) *= column( ref_              , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( max( olhs_, rhs_ ), j );
               column( odres_ , j ) *= column( max( olhs_, rhs_ ), j );
               column( sres_  , j ) *= column( max( olhs_, rhs_ ), j );
               column( osres_ , j ) *= column( max( olhs_, rhs_ ), j );
               column( refres_, j ) *= column( ref_              , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( max( olhs_, orhs_ ), j );
               column( odres_ , j ) *= column( max( olhs_, orhs_ ), j );
               column( sres_  , j ) *= column( max( olhs_, orhs_ ), j );
               column( osres_ , j ) *= column( max( olhs_, orhs_ ), j );
               column( refres_, j ) *= column( ref_               , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Column-wise maximum with multiplication assignment with evaluated matrices
      {
         test_  = "Column-wise maximum with multiplication assignment with evaluated matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( odres_ , j ) *= column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( sres_  , j ) *= column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( osres_ , j ) *= column( max( eval( lhs_ ), eval( rhs_ ) ), j );
               column( refres_, j ) *= column( eval( ref_ )                     , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( odres_ , j ) *= column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( sres_  , j ) *= column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( osres_ , j ) *= column( max( eval( lhs_ ), eval( orhs_ ) ), j );
               column( refres_, j ) *= column( eval( ref_ )                      , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( odres_ , j ) *= column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( sres_  , j ) *= column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( osres_ , j ) *= column( max( eval( olhs_ ), eval( rhs_ ) ), j );
               column( refres_, j ) *= column( eval( ref_ )                      , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<lhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( odres_ , j ) *= column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( sres_  , j ) *= column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( osres_ , j ) *= column( max( eval( olhs_ ), eval( orhs_ ) ), j );
               column( refres_, j ) *= column( eval( ref_ )                       , j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid column index)
      {
         test_  = "Out-of-bounds column construction (invalid column index)";
         error_ = "Setup of out-of-bounds column succeeded";

         try {
            auto c = column( max( lhs_, rhs_ ), lhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = column( max( lhs_, orhs_ ), lhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = column( max( olhs_, rhs_ ), olhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = column( max( olhs_, orhs_ ), olhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
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
/*!\brief Skipping the column-wise dense matrix/dense matrix maximum operation.
//
// \return void
//
// This function is called in case the column-wise matrix/matrix maximum operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testColumnOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the columns-wise dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the columns-wise matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testColumnsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION > 1 )
   {
      if( lhs_.columns() == 0UL )
         return;


      std::vector<size_t> indices( lhs_.columns() );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );


      //=====================================================================================
      // Columns-wise maximum
      //=====================================================================================

      // Columns-wise maximum with the given matrices
      {
         test_  = "Columns-wise maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Columns-wise maximum with evaluated matrices
      {
         test_  = "Columns-wise maximum with evaluated matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Columns-wise maximum with addition assignment
      //=====================================================================================

      // Columns-wise maximum with addition assignment with the given matrices
      {
         test_  = "Columns-wise maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Columns-wise maximum with addition assignment with evaluated matrices
      {
         test_  = "Columns-wise maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Columns-wise maximum with subtraction assignment
      //=====================================================================================

      // Columns-wise maximum with subtraction assignment with the given matrices
      {
         test_  = "Columns-wise maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Columns-wise maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Columns-wise maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Columns-wise maximum with Schur product assignment
      //=====================================================================================

      // Columns-wise maximum with Schur product assignment with the given matrices
      {
         test_  = "Columns-wise maximum with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( max( lhs_, rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( max( lhs_, orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( max( olhs_, rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( max( olhs_, orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( ref_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Columns-wise maximum with Schur product assignment with evaluated matrices
      {
         test_  = "Columns-wise maximum with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( max( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( max( eval( lhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( max( eval( olhs_ ), eval( rhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( max( eval( olhs_ ), eval( orhs_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( eval( ref_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid column index; initializer_list)
      {
         test_  = "Out-of-bounds column selection construction (invalid column index; initializer_list)";
         error_ = "Setup of out-of-bounds column selection succeeded";

         try {
            auto c = columns( max( lhs_, rhs_ ), { lhs_.columns() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( max( lhs_, orhs_ ), { lhs_.columns() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( max( olhs_, rhs_ ), { olhs_.columns() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( max( olhs_, orhs_ ), { olhs_.columns() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
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
            auto c = columns( max( lhs_, rhs_ ), [index=lhs_.columns()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( max( lhs_, orhs_ ), [index=lhs_.columns()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( max( olhs_, rhs_ ), [index=olhs_.columns()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( max( olhs_, orhs_ ), [index=olhs_.columns()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
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
/*!\brief Skipping the columns-wise dense matrix/dense matrix maximum operation.
//
// \return void
//
// This function is called in case the columns-wise matrix/matrix maximum operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testColumnsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the band-wise dense matrix/dense matrix maximum operation.
//
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the band-wise matrix maximum with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// maximum operation or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testBandOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_BAND_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BAND_OPERATION > 1 )
   {
      if( lhs_.rows() == 0UL || lhs_.columns() == 0UL )
         return;


      const ptrdiff_t ibegin( 1UL - lhs_.rows() );
      const ptrdiff_t iend  ( lhs_.columns() );


      //=====================================================================================
      // Band-wise maximum
      //=====================================================================================

      // Band-wise maximum with the given matrices
      {
         test_  = "Band-wise maximum with the given matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( max( lhs_, rhs_ ), i );
               band( odres_ , i ) = band( max( lhs_, rhs_ ), i );
               band( sres_  , i ) = band( max( lhs_, rhs_ ), i );
               band( osres_ , i ) = band( max( lhs_, rhs_ ), i );
               band( refres_, i ) = band( ref_             , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( max( lhs_, orhs_ ), i );
               band( odres_ , i ) = band( max( lhs_, orhs_ ), i );
               band( sres_  , i ) = band( max( lhs_, orhs_ ), i );
               band( osres_ , i ) = band( max( lhs_, orhs_ ), i );
               band( refres_, i ) = band( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( max( olhs_, rhs_ ), i );
               band( odres_ , i ) = band( max( olhs_, rhs_ ), i );
               band( sres_  , i ) = band( max( olhs_, rhs_ ), i );
               band( osres_ , i ) = band( max( olhs_, rhs_ ), i );
               band( refres_, i ) = band( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( max( olhs_, orhs_ ), i );
               band( odres_ , i ) = band( max( olhs_, orhs_ ), i );
               band( sres_  , i ) = band( max( olhs_, orhs_ ), i );
               band( osres_ , i ) = band( max( olhs_, orhs_ ), i );
               band( refres_, i ) = band( ref_               , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Band-wise maximum with evaluated matrices
      {
         test_  = "Band-wise maximum with evaluated matrices";
         error_ = "Failed maximum operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( odres_ , i ) = band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( sres_  , i ) = band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( osres_ , i ) = band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( refres_, i ) = band( eval( ref_ )                     , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( odres_ , i ) = band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( sres_  , i ) = band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( osres_ , i ) = band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( refres_, i ) = band( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( odres_ , i ) = band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( sres_  , i ) = band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( osres_ , i ) = band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( refres_, i ) = band( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( odres_ , i ) = band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( sres_  , i ) = band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( osres_ , i ) = band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( refres_, i ) = band( eval( ref_ )                       , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Band-wise maximum with addition assignment
      //=====================================================================================

      // Band-wise maximum with addition assignment with the given matrices
      {
         test_  = "Band-wise maximum with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( max( lhs_, rhs_ ), i );
               band( odres_ , i ) += band( max( lhs_, rhs_ ), i );
               band( sres_  , i ) += band( max( lhs_, rhs_ ), i );
               band( osres_ , i ) += band( max( lhs_, rhs_ ), i );
               band( refres_, i ) += band( ref_             , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( max( lhs_, orhs_ ), i );
               band( odres_ , i ) += band( max( lhs_, orhs_ ), i );
               band( sres_  , i ) += band( max( lhs_, orhs_ ), i );
               band( osres_ , i ) += band( max( lhs_, orhs_ ), i );
               band( refres_, i ) += band( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( max( olhs_, rhs_ ), i );
               band( odres_ , i ) += band( max( olhs_, rhs_ ), i );
               band( sres_  , i ) += band( max( olhs_, rhs_ ), i );
               band( osres_ , i ) += band( max( olhs_, rhs_ ), i );
               band( refres_, i ) += band( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( max( olhs_, orhs_ ), i );
               band( odres_ , i ) += band( max( olhs_, orhs_ ), i );
               band( sres_  , i ) += band( max( olhs_, orhs_ ), i );
               band( osres_ , i ) += band( max( olhs_, orhs_ ), i );
               band( refres_, i ) += band( ref_               , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Band-wise maximum with addition assignment with evaluated matrices
      {
         test_  = "Band-wise maximum with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( odres_ , i ) += band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( sres_  , i ) += band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( osres_ , i ) += band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( refres_, i ) += band( eval( ref_ )                     , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( odres_ , i ) += band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( sres_  , i ) += band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( osres_ , i ) += band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( refres_, i ) += band( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( odres_ , i ) += band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( sres_  , i ) += band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( osres_ , i ) += band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( refres_, i ) += band( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( odres_ , i ) += band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( sres_  , i ) += band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( osres_ , i ) += band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( refres_, i ) += band( eval( ref_ )                       , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Band-wise maximum with subtraction assignment
      //=====================================================================================

      // Band-wise maximum with subtraction assignment with the given matrices
      {
         test_  = "Band-wise maximum with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( max( lhs_, rhs_ ), i );
               band( odres_ , i ) -= band( max( lhs_, rhs_ ), i );
               band( sres_  , i ) -= band( max( lhs_, rhs_ ), i );
               band( osres_ , i ) -= band( max( lhs_, rhs_ ), i );
               band( refres_, i ) -= band( ref_             , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( max( lhs_, orhs_ ), i );
               band( odres_ , i ) -= band( max( lhs_, orhs_ ), i );
               band( sres_  , i ) -= band( max( lhs_, orhs_ ), i );
               band( osres_ , i ) -= band( max( lhs_, orhs_ ), i );
               band( refres_, i ) -= band( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( max( olhs_, rhs_ ), i );
               band( odres_ , i ) -= band( max( olhs_, rhs_ ), i );
               band( sres_  , i ) -= band( max( olhs_, rhs_ ), i );
               band( osres_ , i ) -= band( max( olhs_, rhs_ ), i );
               band( refres_, i ) -= band( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( max( olhs_, orhs_ ), i );
               band( odres_ , i ) -= band( max( olhs_, orhs_ ), i );
               band( sres_  , i ) -= band( max( olhs_, orhs_ ), i );
               band( osres_ , i ) -= band( max( olhs_, orhs_ ), i );
               band( refres_, i ) -= band( ref_               , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Band-wise maximum with subtraction assignment with evaluated matrices
      {
         test_  = "Band-wise maximum with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( odres_ , i ) -= band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( sres_  , i ) -= band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( osres_ , i ) -= band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( refres_, i ) -= band( eval( ref_ )                     , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( odres_ , i ) -= band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( sres_  , i ) -= band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( osres_ , i ) -= band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( refres_, i ) -= band( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( odres_ , i ) -= band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( sres_  , i ) -= band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( osres_ , i ) -= band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( refres_, i ) -= band( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( odres_ , i ) -= band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( sres_  , i ) -= band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( osres_ , i ) -= band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( refres_, i ) -= band( eval( ref_ )                       , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Band-wise maximum with multiplication assignment
      //=====================================================================================

      // Band-wise maximum with multiplication assignment with the given matrices
      {
         test_  = "Band-wise maximum with multiplication assignment with the given matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( max( lhs_, rhs_ ), i );
               band( odres_ , i ) *= band( max( lhs_, rhs_ ), i );
               band( sres_  , i ) *= band( max( lhs_, rhs_ ), i );
               band( osres_ , i ) *= band( max( lhs_, rhs_ ), i );
               band( refres_, i ) *= band( ref_             , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( max( lhs_, orhs_ ), i );
               band( odres_ , i ) *= band( max( lhs_, orhs_ ), i );
               band( sres_  , i ) *= band( max( lhs_, orhs_ ), i );
               band( osres_ , i ) *= band( max( lhs_, orhs_ ), i );
               band( refres_, i ) *= band( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( max( olhs_, rhs_ ), i );
               band( odres_ , i ) *= band( max( olhs_, rhs_ ), i );
               band( sres_  , i ) *= band( max( olhs_, rhs_ ), i );
               band( osres_ , i ) *= band( max( olhs_, rhs_ ), i );
               band( refres_, i ) *= band( ref_              , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( max( olhs_, orhs_ ), i );
               band( odres_ , i ) *= band( max( olhs_, orhs_ ), i );
               band( sres_  , i ) *= band( max( olhs_, orhs_ ), i );
               band( osres_ , i ) *= band( max( olhs_, orhs_ ), i );
               band( refres_, i ) *= band( ref_               , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Band-wise maximum with multiplication assignment with evaluated matrices
      {
         test_  = "Band-wise maximum with multiplication assignment with evaluated matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( odres_ , i ) *= band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( sres_  , i ) *= band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( osres_ , i ) *= band( max( eval( lhs_ ), eval( rhs_ ) ), i );
               band( refres_, i ) *= band( eval( ref_ )                     , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( odres_ , i ) *= band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( sres_  , i ) *= band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( osres_ , i ) *= band( max( eval( lhs_ ), eval( orhs_ ) ), i );
               band( refres_, i ) *= band( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( odres_ , i ) *= band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( sres_  , i ) *= band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( osres_ , i ) *= band( max( eval( olhs_ ), eval( rhs_ ) ), i );
               band( refres_, i ) *= band( eval( ref_ )                      , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( odres_ , i ) *= band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( sres_  , i ) *= band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( osres_ , i ) *= band( max( eval( olhs_ ), eval( orhs_ ) ), i );
               band( refres_, i ) *= band( eval( ref_ )                       , i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid lower band index)
      {
         test_  = "Out-of-bounds band construction (invalid lower band index)";
         error_ = "Setup of out-of-bounds band succeeded";

         try {
            auto b = band( max( lhs_, rhs_ ), -lhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( max( lhs_, orhs_ ), -lhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( max( olhs_, rhs_ ), -olhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( max( olhs_, orhs_ ), -olhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
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
            auto b = band( max( lhs_, rhs_ ), lhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( max( lhs_, orhs_ ), lhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( max( olhs_, rhs_ ), olhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( max( olhs_, orhs_ ), olhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side dense matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side dense matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
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
/*!\brief Skipping the band-wise dense matrix/dense matrix maximum operation.
//
// \return void
//
// This function is called in case the band-wise matrix/matrix maximum operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::testBandOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the customized dense matrix/dense matrix maximum operation.
//
// \param op The custom operation to be tested.
// \param name The human-readable name of the operation.
// \return void
// \exception std::runtime_error Maximum error detected.
//
// This function tests the matrix maximum with plain assignment, addition assignment, and
// subtraction assignment in combination with a custom operation. In case any error resulting
// from the maximum operation or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
template< typename OP >   // Type of the custom operation
void OperationTest<MT1,MT2>::testCustomOperation( OP op, const std::string& name )
{
   //=====================================================================================
   // Customized maximum
   //=====================================================================================

   // Customized maximum with the given matrices
   {
      test_  = "Customized maximum with the given matrices (" + name + ")";
      error_ = "Failed maximum operation";

      try {
         initResults();
         dres_   = op( max( lhs_, rhs_ ) );
         odres_  = op( max( lhs_, rhs_ ) );
         sres_   = op( max( lhs_, rhs_ ) );
         osres_  = op( max( lhs_, rhs_ ) );
         refres_ = op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   = op( max( lhs_, orhs_ ) );
         odres_  = op( max( lhs_, orhs_ ) );
         sres_   = op( max( lhs_, orhs_ ) );
         osres_  = op( max( lhs_, orhs_ ) );
         refres_ = op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   = op( max( olhs_, rhs_ ) );
         odres_  = op( max( olhs_, rhs_ ) );
         sres_   = op( max( olhs_, rhs_ ) );
         osres_  = op( max( olhs_, rhs_ ) );
         refres_ = op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   = op( max( olhs_, orhs_ ) );
         odres_  = op( max( olhs_, orhs_ ) );
         sres_   = op( max( olhs_, orhs_ ) );
         osres_  = op( max( olhs_, orhs_ ) );
         refres_ = op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }

   // Customized maximum with evaluated matrices
   {
      test_  = "Customized maximum with evaluated matrices (" + name + ")";
      error_ = "Failed maximum operation";

      try {
         initResults();
         dres_   = op( max( eval( lhs_ ), eval( rhs_ ) ) );
         odres_  = op( max( eval( lhs_ ), eval( rhs_ ) ) );
         sres_   = op( max( eval( lhs_ ), eval( rhs_ ) ) );
         osres_  = op( max( eval( lhs_ ), eval( rhs_ ) ) );
         refres_ = op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   = op( max( eval( lhs_ ), eval( orhs_ ) ) );
         odres_  = op( max( eval( lhs_ ), eval( orhs_ ) ) );
         sres_   = op( max( eval( lhs_ ), eval( orhs_ ) ) );
         osres_  = op( max( eval( lhs_ ), eval( orhs_ ) ) );
         refres_ = op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   = op( max( eval( olhs_ ), eval( rhs_ ) ) );
         odres_  = op( max( eval( olhs_ ), eval( rhs_ ) ) );
         sres_   = op( max( eval( olhs_ ), eval( rhs_ ) ) );
         osres_  = op( max( eval( olhs_ ), eval( rhs_ ) ) );
         refres_ = op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   = op( max( eval( olhs_ ), eval( orhs_ ) ) );
         odres_  = op( max( eval( olhs_ ), eval( orhs_ ) ) );
         sres_   = op( max( eval( olhs_ ), eval( orhs_ ) ) );
         osres_  = op( max( eval( olhs_ ), eval( orhs_ ) ) );
         refres_ = op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }


   //=====================================================================================
   // Customized maximum with addition assignment
   //=====================================================================================

   // Customized maximum with addition assignment with the given matrices
   {
      test_  = "Customized maximum with addition assignment with the given matrices (" + name + ")";
      error_ = "Failed addition assignment operation";

      try {
         initResults();
         dres_   += op( max( lhs_, rhs_ ) );
         odres_  += op( max( lhs_, rhs_ ) );
         sres_   += op( max( lhs_, rhs_ ) );
         osres_  += op( max( lhs_, rhs_ ) );
         refres_ += op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   += op( max( lhs_, orhs_ ) );
         odres_  += op( max( lhs_, orhs_ ) );
         sres_   += op( max( lhs_, orhs_ ) );
         osres_  += op( max( lhs_, orhs_ ) );
         refres_ += op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   += op( max( olhs_, rhs_ ) );
         odres_  += op( max( olhs_, rhs_ ) );
         sres_   += op( max( olhs_, rhs_ ) );
         osres_  += op( max( olhs_, rhs_ ) );
         refres_ += op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   += op( max( olhs_, orhs_ ) );
         odres_  += op( max( olhs_, orhs_ ) );
         sres_   += op( max( olhs_, orhs_ ) );
         osres_  += op( max( olhs_, orhs_ ) );
         refres_ += op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }

   // Customized maximum with addition assignment with evaluated matrices
   {
      test_  = "Customized maximum with addition assignment with evaluated matrices (" + name + ")";
      error_ = "Failed addition assignment operation";

      try {
         initResults();
         dres_   += op( max( eval( lhs_ ), eval( rhs_ ) ) );
         odres_  += op( max( eval( lhs_ ), eval( rhs_ ) ) );
         sres_   += op( max( eval( lhs_ ), eval( rhs_ ) ) );
         osres_  += op( max( eval( lhs_ ), eval( rhs_ ) ) );
         refres_ += op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   += op( max( eval( lhs_ ), eval( orhs_ ) ) );
         odres_  += op( max( eval( lhs_ ), eval( orhs_ ) ) );
         sres_   += op( max( eval( lhs_ ), eval( orhs_ ) ) );
         osres_  += op( max( eval( lhs_ ), eval( orhs_ ) ) );
         refres_ += op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   += op( max( eval( olhs_ ), eval( rhs_ ) ) );
         odres_  += op( max( eval( olhs_ ), eval( rhs_ ) ) );
         sres_   += op( max( eval( olhs_ ), eval( rhs_ ) ) );
         osres_  += op( max( eval( olhs_ ), eval( rhs_ ) ) );
         refres_ += op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   += op( max( eval( olhs_ ), eval( orhs_ ) ) );
         odres_  += op( max( eval( olhs_ ), eval( orhs_ ) ) );
         sres_   += op( max( eval( olhs_ ), eval( orhs_ ) ) );
         osres_  += op( max( eval( olhs_ ), eval( orhs_ ) ) );
         refres_ += op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }


   //=====================================================================================
   // Customized maximum with subtraction assignment
   //=====================================================================================

   // Customized maximum with subtraction assignment with the given matrices
   {
      test_  = "Customized maximum with subtraction assignment with the given matrices (" + name + ")";
      error_ = "Failed subtraction assignment operation";

      try {
         initResults();
         dres_   -= op( max( lhs_, rhs_ ) );
         odres_  -= op( max( lhs_, rhs_ ) );
         sres_   -= op( max( lhs_, rhs_ ) );
         osres_  -= op( max( lhs_, rhs_ ) );
         refres_ -= op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   -= op( max( lhs_, orhs_ ) );
         odres_  -= op( max( lhs_, orhs_ ) );
         sres_   -= op( max( lhs_, orhs_ ) );
         osres_  -= op( max( lhs_, orhs_ ) );
         refres_ -= op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   -= op( max( olhs_, rhs_ ) );
         odres_  -= op( max( olhs_, rhs_ ) );
         sres_   -= op( max( olhs_, rhs_ ) );
         osres_  -= op( max( olhs_, rhs_ ) );
         refres_ -= op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   -= op( max( olhs_, orhs_ ) );
         odres_  -= op( max( olhs_, orhs_ ) );
         sres_   -= op( max( olhs_, orhs_ ) );
         osres_  -= op( max( olhs_, orhs_ ) );
         refres_ -= op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }

   // Customized maximum with subtraction assignment with evaluated matrices
   {
      test_  = "Customized maximum with subtraction assignment with evaluated matrices (" + name + ")";
      error_ = "Failed subtraction assignment operation";

      try {
         initResults();
         dres_   -= op( max( eval( lhs_ ), eval( rhs_ ) ) );
         odres_  -= op( max( eval( lhs_ ), eval( rhs_ ) ) );
         sres_   -= op( max( eval( lhs_ ), eval( rhs_ ) ) );
         osres_  -= op( max( eval( lhs_ ), eval( rhs_ ) ) );
         refres_ -= op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   -= op( max( eval( lhs_ ), eval( orhs_ ) ) );
         odres_  -= op( max( eval( lhs_ ), eval( orhs_ ) ) );
         sres_   -= op( max( eval( lhs_ ), eval( orhs_ ) ) );
         osres_  -= op( max( eval( lhs_ ), eval( orhs_ ) ) );
         refres_ -= op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   -= op( max( eval( olhs_ ), eval( rhs_ ) ) );
         odres_  -= op( max( eval( olhs_ ), eval( rhs_ ) ) );
         sres_   -= op( max( eval( olhs_ ), eval( rhs_ ) ) );
         osres_  -= op( max( eval( olhs_ ), eval( rhs_ ) ) );
         refres_ -= op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   -= op( max( eval( olhs_ ), eval( orhs_ ) ) );
         odres_  -= op( max( eval( olhs_ ), eval( orhs_ ) ) );
         sres_   -= op( max( eval( olhs_ ), eval( orhs_ ) ) );
         osres_  -= op( max( eval( olhs_ ), eval( orhs_ ) ) );
         refres_ -= op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }


   //=====================================================================================
   // Customized maximum with Schur product assignment
   //=====================================================================================

   // Customized maximum with Schur product assignment with the given matrices
   {
      test_  = "Customized maximum with Schur product assignment with the given matrices (" + name + ")";
      error_ = "Failed Schur product assignment operation";

      try {
         initResults();
         dres_   %= op( max( lhs_, rhs_ ) );
         odres_  %= op( max( lhs_, rhs_ ) );
         sres_   %= op( max( lhs_, rhs_ ) );
         osres_  %= op( max( lhs_, rhs_ ) );
         refres_ %= op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   %= op( max( lhs_, orhs_ ) );
         odres_  %= op( max( lhs_, orhs_ ) );
         sres_   %= op( max( lhs_, orhs_ ) );
         osres_  %= op( max( lhs_, orhs_ ) );
         refres_ %= op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   %= op( max( olhs_, rhs_ ) );
         odres_  %= op( max( olhs_, rhs_ ) );
         sres_   %= op( max( olhs_, rhs_ ) );
         osres_  %= op( max( olhs_, rhs_ ) );
         refres_ %= op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   %= op( max( olhs_, orhs_ ) );
         odres_  %= op( max( olhs_, orhs_ ) );
         sres_   %= op( max( olhs_, orhs_ ) );
         osres_  %= op( max( olhs_, orhs_ ) );
         refres_ %= op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }

   // Customized maximum with Schur product assignment with evaluated matrices
   {
      test_  = "Customized maximum with Schur product assignment with evaluated matrices (" + name + ")";
      error_ = "Failed Schur product assignment operation";

      try {
         initResults();
         dres_   %= op( max( eval( lhs_ ), eval( rhs_ ) ) );
         odres_  %= op( max( eval( lhs_ ), eval( rhs_ ) ) );
         sres_   %= op( max( eval( lhs_ ), eval( rhs_ ) ) );
         osres_  %= op( max( eval( lhs_ ), eval( rhs_ ) ) );
         refres_ %= op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   %= op( max( eval( lhs_ ), eval( orhs_ ) ) );
         odres_  %= op( max( eval( lhs_ ), eval( orhs_ ) ) );
         sres_   %= op( max( eval( lhs_ ), eval( orhs_ ) ) );
         osres_  %= op( max( eval( lhs_ ), eval( orhs_ ) ) );
         refres_ %= op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   %= op( max( eval( olhs_ ), eval( rhs_ ) ) );
         odres_  %= op( max( eval( olhs_ ), eval( rhs_ ) ) );
         sres_   %= op( max( eval( olhs_ ), eval( rhs_ ) ) );
         osres_  %= op( max( eval( olhs_ ), eval( rhs_ ) ) );
         refres_ %= op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   %= op( max( eval( olhs_ ), eval( orhs_ ) ) );
         odres_  %= op( max( eval( olhs_ ), eval( orhs_ ) ) );
         sres_   %= op( max( eval( olhs_ ), eval( orhs_ ) ) );
         osres_  %= op( max( eval( olhs_ ), eval( orhs_ ) ) );
         refres_ %= op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
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
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
template< typename LT     // Type of the left-hand side operand
        , typename RT >   // Type of the right-hand side operand
void OperationTest<MT1,MT2>::checkResults()
{
   using blaze::IsRowMajorMatrix;

   if( !isEqual( dres_, refres_ ) || !isEqual( odres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( RT ).name() << "\n"
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
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( RT ).name() << "\n"
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
// results. The two template arguments \a LT and \a RT indicate the types of the left-hand
// side and right-hand side operands used for the computations.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
template< typename LT     // Type of the left-hand side operand
        , typename RT >   // Type of the right-hand side operand
void OperationTest<MT1,MT2>::checkTransposeResults()
{
   using blaze::IsRowMajorMatrix;

   if( !isEqual( tdres_, refres_ ) || !isEqual( todres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( RT ).name() << "\n"
          << "   Transpose result:\n" << tdres_ << "\n"
          << "   Transpose result with opposite storage order:\n" << todres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( tsres_, refres_ ) || !isEqual( tosres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
          << "     " << typeid( RT ).name() << "\n"
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
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::checkExceptionMessage( const std::exception& ex, const std::string& message )
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
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::initResults()
{
   const blaze::UnderlyingBuiltin_t<DRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<DRE> max( randmax );

   resize( dres_, rows( lhs_ ), columns( lhs_ ) );
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
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void OperationTest<MT1,MT2>::initTransposeResults()
{
   const blaze::UnderlyingBuiltin_t<TDRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<TDRE> max( randmax );

   resize( tdres_, columns( lhs_ ), rows( lhs_ ) );
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
// test. The two template arguments \a LT and \a RT indicate the types of the left-hand side and
// right-hand side operands used for the computations.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
template< typename LT     // Type of the left-hand side operand
        , typename RT >   // Type of the right-hand side operand
void OperationTest<MT1,MT2>::convertException( const std::exception& ex )
{
   using blaze::IsRowMajorMatrix;

   std::ostringstream oss;
   oss << " Test : " << test_ << "\n"
       << " Error: " << error_ << "\n"
       << " Details:\n"
       << "   Random seed = " << blaze::getSeed() << "\n"
       << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
       << "     " << typeid( LT ).name() << "\n"
       << "   Right-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " dense matrix type:\n"
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
/*!\brief Testing the matrix maximum operation between two specific matrix types.
//
// \param creator1 The creator for the left-hand side matrix.
// \param creator2 The creator for the right-hand side matrix.
// \return void
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
void runTest( const Creator<MT1>& creator1, const Creator<MT2>& creator2 )
{
#if BLAZETEST_MATHTEST_TEST_MAXIMUM
   if( BLAZETEST_MATHTEST_TEST_MAXIMUM > 1 )
   {
      for( size_t rep=0UL; rep<BLAZETEST_REPETITIONS; ++rep ) {
         OperationTest<MT1,MT2>( creator1, creator2 );
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
/*!\brief Macro for the definition of a dense matrix/dense matrix maximum test case.
*/
#define DEFINE_DMATDMATMAX_OPERATION_TEST( MT1, MT2 ) \
   extern template class blazetest::mathtest::operations::dmatdmatmax::OperationTest<MT1,MT2>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a dense matrix/dense matrix maximum test case.
*/
#define RUN_DMATDMATMAX_OPERATION_TEST( C1, C2 ) \
   blazetest::mathtest::operations::dmatdmatmax::runTest( C1, C2 )
/*! \endcond */
//*************************************************************************************************

} // namespace dmatdmatmax

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
