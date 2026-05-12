//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/smatsmatmult/OperationTest.h
//  \brief Header file for the sparse matrix/sparse matrix multiplication operation test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_SMATSMATMULT_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_SMATSMATMULT_OPERATIONTEST_H_


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
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingScalar.h>
#include <blaze/math/Views.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/constraints/DerivedFrom.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Nor.h>
#include <blaze/util/mpl/Or.h>
#include <blaze/util/Random.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveCVRef.h>
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

namespace smatsmatmult {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the sparse matrix/sparse matrix multiplication operation test.
//
// This class template represents one particular matrix multiplication test between two matrices
// of a particular type. The two template arguments \a MT1 and \a MT2 represent the types of the
// left-hand side and right-hand side matrix, respectively.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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

   //! Sparse result type
   using SRE = blaze::MultTrait_t<MT1,MT2>;

   using SET   = blaze::ElementType_t<SRE>;     //!< Element type of the sparse result
   using OSRE  = blaze::OppositeType_t<SRE>;    //!< Sparse result type with opposite storage order
   using TSRE  = blaze::TransposeType_t<SRE>;   //!< Transpose sparse result type
   using TOSRE = blaze::TransposeType_t<OSRE>;  //!< Transpose sparse result type with opposite storage order

   //! Dense result type
   using DRE = MatchAdaptor_t< SRE, blaze::DynamicMatrix<SET,false> >;

   using DET   = blaze::ElementType_t<DRE>;     //!< Element type of the dense result
   using ODRE  = blaze::OppositeType_t<DRE>;    //!< Dense result type with opposite storage order
   using TDRE  = blaze::TransposeType_t<DRE>;   //!< Transpose dense result type
   using TODRE = blaze::TransposeType_t<ODRE>;  //!< Transpose dense result type with opposite storage order

   using RT1 = blaze::DynamicMatrix<ET1,false>;  //!< Reference type 1
   using RT2 = blaze::DynamicMatrix<ET2,false>;  //!< Reference type 2

   //! Reference result type
   using RRE = MatchSymmetry_t< DRE, blaze::MultTrait_t<RT1,RT2> >;

   //! Type of the matrix/matrix multiplication expression
   using MatMatMultExprType =
      blaze::RemoveCVRef_t< decltype( std::declval<MT1>() * std::declval<MT2>() ) >;

   //! Type of the matrix/transpose matrix multiplication expression
   using MatTMatMultExprType =
      blaze::RemoveCVRef_t< decltype( std::declval<MT1>() * std::declval<OMT2>() ) >;

   //! Type of the transpose matrix/matrix multiplication expression
   using TMatMatMultExprType =
      blaze::RemoveCVRef_t< decltype( std::declval<OMT1>() * std::declval<MT2>() ) >;

   //! Type of the transpose matrix/transpose matrix multiplication expression
   using TMatTMatMultExprType =
      blaze::RemoveCVRef_t< decltype( std::declval<OMT1>() * std::declval<OMT2>() ) >;
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
                          void testSubmatrixOperation();
                          void testRowOperation      ();
                          void testRowsOperation     ( blaze::TrueType  );
                          void testRowsOperation     ( blaze::FalseType );
                          void testColumnOperation   ();
                          void testColumnsOperation  ( blaze::TrueType  );
                          void testColumnsOperation  ( blaze::FalseType );
                          void testBandOperation     ();

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
   MT1   lhs_;     //!< The left-hand side sparse matrix.
   MT2   rhs_;     //!< The right-hand side sparse matrix.
   OMT1  olhs_;    //!< The left-hand side sparse matrix with opposite storage order.
   OMT2  orhs_;    //!< The right-hand side sparse matrix with opposite storage order.
   DRE   dres_;    //!< The dense result matrix.
   SRE   sres_;    //!< The sparse result matrix.
   ODRE  odres_;   //!< The dense result matrix with opposite storage order.
   OSRE  osres_;   //!< The sparse result matrix with opposite storage order.
   TDRE  tdres_;   //!< The transpose dense result matrix.
   TSRE  tsres_;   //!< The transpose sparse result matrix.
   TODRE todres_;  //!< The transpose dense result matrix with opposite storage order.
   TOSRE tosres_;  //!< The transpose sparse result matrix with opposite storage order.
   RT1   reflhs_;  //!< The reference left-hand side matrix.
   RT2   refrhs_;  //!< The reference right-hand side matrix.
   RRE   refres_;  //!< The reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT1   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT2   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OMT1  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OMT2  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TMT1  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TMT2  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TOMT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TOMT2 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( RT1   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( RT2   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( RRE   );
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
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( RT1   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( RT2   );
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

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( MatMatMultExprType, blaze::ResultType_t<MatMatMultExprType>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MatMatMultExprType, blaze::OppositeType_t<MatMatMultExprType>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MatMatMultExprType, blaze::TransposeType_t<MatMatMultExprType> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( MatTMatMultExprType, blaze::ResultType_t<MatTMatMultExprType>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MatTMatMultExprType, blaze::OppositeType_t<MatTMatMultExprType>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MatTMatMultExprType, blaze::TransposeType_t<MatTMatMultExprType> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( TMatMatMultExprType, blaze::ResultType_t<TMatMatMultExprType>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TMatMatMultExprType, blaze::OppositeType_t<TMatMatMultExprType>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TMatMatMultExprType, blaze::TransposeType_t<TMatMatMultExprType> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( TMatTMatMultExprType, blaze::ResultType_t<TMatTMatMultExprType>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TMatTMatMultExprType, blaze::OppositeType_t<TMatTMatMultExprType>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TMatTMatMultExprType, blaze::TransposeType_t<TMatTMatMultExprType> );

   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( MatMatMultExprType  , blaze::BaseType_t<MatMatMultExprType  > );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( MatTMatMultExprType , blaze::BaseType_t<MatTMatMultExprType > );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TMatMatMultExprType , blaze::BaseType_t<TMatMatMultExprType > );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TMatTMatMultExprType, blaze::BaseType_t<TMatTMatMultExprType> );
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
/*!\brief Constructor for the sparse matrix/sparse matrix multiplication operation test.
//
// \param creator1 The creator for the left-hand side sparse matrix of the matrix multiplication.
// \param creator2 The creator for the right-hand side sparse matrix of the matrix multiplication.
// \exception std::runtime_error Operation error detected.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
OperationTest<MT1,MT2>::OperationTest( const Creator<MT1>& creator1, const Creator<MT2>& creator2 )
   : lhs_( creator1() )  // The left-hand side sparse matrix
   , rhs_( creator2() )  // The right-hand side sparse matrix
   , olhs_( lhs_ )       // The left-hand side sparse matrix with opposite storage order
   , orhs_( rhs_ )       // The right-hand side sparse matrix with opposite storage order
   , dres_()             // The dense result matrix
   , sres_()             // The sparse result matrix
   , odres_()            // The dense result matrix with opposite storage order
   , osres_()            // The sparse result matrix with opposite storage order
   , tdres_()            // The transpose dense result matrix
   , tsres_()            // The transpose sparse result matrix
   , todres_()           // The transpose dense result matrix with opposite storage order
   , tosres_()           // The transpose sparse result matrix with opposite storage order
   , reflhs_( lhs_ )     // The reference left-hand side matrix
   , refrhs_( rhs_ )     // The reference right-hand side matrix
   , refres_()           // The reference result
   , test_()             // Label of the currently performed test
   , error_()            // Description of the current error type
{
   using namespace blaze;

   using Scalar = UnderlyingScalar_t<SET>;

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
   testDeclSymOperation( IsSame< LowType_t<ET1,ET2>, ET1 >() );
   testDeclHermOperation( IsSame< LowType_t<ET1,ET2>, ET1 >() );
   testDeclLowOperation( Or_t< IsSquare<SRE>, IsResizable<SRE> >() );
   testDeclUppOperation( Or_t< IsSquare<SRE>, IsResizable<SRE> >() );
   testDeclDiagOperation( Or_t< IsSquare<SRE>, IsResizable<SRE> >() );
   testSubmatrixOperation();
   testRowOperation();
   testRowsOperation( Nor_t< IsSymmetric<DRE>, IsHermitian<DRE> >() );
   testColumnOperation();
   testColumnsOperation( Nor_t< IsSymmetric<DRE>, IsHermitian<DRE> >() );
   testBandOperation();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Tests on the initial status of the matrices.
//
// \return void
// \exception std::runtime_error Initialization error detected.
//
// This function runs tests on the initial status of the matrices. In case any initialization
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testInitialStatus()
{
   //=====================================================================================
   // Performing initial tests with the row-major types
   //=====================================================================================

   // Checking the number of rows of the left-hand side operand
   if( lhs_.rows() != reflhs_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side row-major sparse operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Detected number of rows = " << lhs_.rows() << "\n"
          << "   Expected number of rows = " << reflhs_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the left-hand side operand
   if( lhs_.columns() != reflhs_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side row-major sparse operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Detected number of columns = " << lhs_.columns() << "\n"
          << "   Expected number of columns = " << reflhs_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of rows of the right-hand side operand
   if( rhs_.rows() != refrhs_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side row-major sparse operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n"
          << "   Detected number of rows = " << rhs_.rows() << "\n"
          << "   Expected number of rows = " << refrhs_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the right-hand side operand
   if( rhs_.columns() != refrhs_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side row-major sparse operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n"
          << "   Detected number of columns = " << rhs_.columns() << "\n"
          << "   Expected number of columns = " << refrhs_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the left-hand side operand
   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of left-hand side row-major sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Current initialization:\n" << lhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the right-hand side operand
   if( !isEqual( rhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of right-hand side row-major sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing initial tests with the column-major types
   //=====================================================================================

   // Checking the number of rows of the left-hand side operand
   if( olhs_.rows() != reflhs_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side column-major dense operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Detected number of rows = " << olhs_.rows() << "\n"
          << "   Expected number of rows = " << reflhs_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the left-hand side operand
   if( olhs_.columns() != reflhs_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of left-hand side column-major dense operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Detected number of columns = " << olhs_.columns() << "\n"
          << "   Expected number of columns = " << reflhs_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of rows of the right-hand side operand
   if( orhs_.rows() != refrhs_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side column-major sparse operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n"
          << "   Detected number of rows = " << orhs_.rows() << "\n"
          << "   Expected number of rows = " << refrhs_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the right-hand side operand
   if( orhs_.columns() != refrhs_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of right-hand side column-major sparse operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n"
          << "   Detected number of columns = " << orhs_.columns() << "\n"
          << "   Expected number of columns = " << refrhs_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the left-hand side operand
   if( !isEqual( olhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of left-hand side column-major dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Current initialization:\n" << olhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the right-hand side operand
   if( !isEqual( orhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of right-hand side column-major sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n"
          << "   Current initialization:\n" << orhs_ << "\n"
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
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testAssignment()
{
   //=====================================================================================
   // // Performing an assignment with the row-major types
   //=====================================================================================

   try {
      lhs_ = reflhs_;
      rhs_ = refrhs_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the row-major types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side row-major sparse matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Right-hand side row-major sparse matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( lhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of left-hand side row-major sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Current initialization:\n" << lhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( rhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of right-hand side row-major sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Row-major sparse matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n"
          << "   Current initialization:\n" << rhs_ << "\n"
          << "   Expected initialization:\n" << refrhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing an assignment with the column-major types
   //=====================================================================================

   try {
      olhs_ = reflhs_;
      orhs_ = refrhs_;
   }
   catch( std::exception& ex ) {
      std::ostringstream oss;
      oss << " Test: Assignment with the column-major types\n"
          << " Error: Failed assignment\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side column-major dense matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Right-hand side column-major sparse matrix type:\n"
          << "     "  << typeid( OMT2 ).name() << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( olhs_, reflhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of left-hand side column-major dense operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Current initialization:\n" << olhs_ << "\n"
          << "   Expected initialization:\n" << reflhs_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( orhs_, refrhs_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of right-hand side column-major sparse operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Column-major sparse matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n"
          << "   Current initialization:\n" << orhs_ << "\n"
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
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testEvaluation()
{
   using blaze::IsRowMajorMatrix;


   //=====================================================================================
   // Testing the evaluation with two row-major matrices
   //=====================================================================================

   {
      const auto res   ( evaluate( lhs_    * rhs_    ) );
      const auto refres( evaluate( reflhs_ * refrhs_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<MT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<MT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
         oss << " Test: Evaluation with evaluated matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<MT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<MT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
      const auto res   ( evaluate( lhs_    * orhs_   ) );
      const auto refres( evaluate( reflhs_ * refrhs_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<MT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<OMT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
      const auto res   ( evaluate( eval( lhs_ )    * eval( orhs_ )   ) );
      const auto refres( evaluate( eval( reflhs_ ) * eval( refrhs_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<MT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
             << "     " << typeid( lhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<OMT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
      const auto res   ( evaluate( olhs_   * rhs_    ) );
      const auto refres( evaluate( reflhs_ * refrhs_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<OMT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
             << "     " << typeid( olhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<MT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<OMT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
             << "     " << typeid( olhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<MT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
      const auto res   ( evaluate( olhs_   * orhs_   ) );
      const auto refres( evaluate( reflhs_ * refrhs_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<OMT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
             << "     " << typeid( olhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<OMT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
      const auto res   ( evaluate( eval( olhs_ )   * eval( orhs_ )   ) );
      const auto refres( evaluate( eval( reflhs_ ) * eval( refrhs_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrices\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side " << ( IsRowMajorMatrix<OMT1>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
             << "     " << typeid( olhs_ ).name() << "\n"
             << "   Right-hand side " << ( IsRowMajorMatrix<OMT2>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testElementAccess()
{
   using blaze::equal;


   //=====================================================================================
   // Testing the element access with two row-major matrices
   //=====================================================================================

   if( lhs_.rows() > 0UL && rhs_.columns() > 0UL )
   {
      const size_t m( lhs_.rows()    - 1UL );
      const size_t n( rhs_.columns() - 1UL );

      if( !equal( ( lhs_ * rhs_ )(m,n), ( reflhs_ * refrhs_ )(m,n) ) ||
          !equal( ( lhs_ * rhs_ ).at(m,n), ( reflhs_ * refrhs_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( lhs_ * eval( rhs_ ) )(m,n), ( reflhs_ * eval( refrhs_ ) )(m,n) ) ||
          !equal( ( lhs_ * eval( rhs_ ) ).at(m,n), ( reflhs_ * eval( refrhs_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * rhs_ )(m,n), ( eval( reflhs_ ) * refrhs_ )(m,n) ) ||
          !equal( ( eval( lhs_ ) * rhs_ ).at(m,n), ( eval( reflhs_ ) * refrhs_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * eval( rhs_ ) )(m,n), ( eval( reflhs_ ) * eval( refrhs_ ) )(m,n) ) ||
          !equal( ( eval( lhs_ ) * eval( rhs_ ) ).at(m,n), ( eval( reflhs_ ) * eval( refrhs_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      ( lhs_ * rhs_ ).at( 0UL, rhs_.columns() );

      std::ostringstream oss;
      oss << " Test : Checked element access of multiplication expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side row-major sparse matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Right-hand side row-major sparse matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      ( lhs_ * rhs_ ).at( lhs_.rows(), 0UL );

      std::ostringstream oss;
      oss << " Test : Checked element access of multiplication expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side row-major sparse matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Right-hand side row-major sparse matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}


   //=====================================================================================
   // Testing the element access with a row-major matrix and a column-major matrix
   //=====================================================================================

   if( lhs_.rows() > 0UL && orhs_.columns() > 0UL )
   {
      const size_t m( lhs_.rows()     - 1UL );
      const size_t n( orhs_.columns() - 1UL );

      if( !equal( ( lhs_ * orhs_ )(m,n), ( reflhs_ * refrhs_ )(m,n) ) ||
          !equal( ( lhs_ * orhs_ ).at(m,n), ( reflhs_ * refrhs_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( lhs_ * eval( orhs_ ) )(m,n), ( reflhs_ * eval( refrhs_ ) )(m,n) ) ||
          !equal( ( lhs_ * eval( orhs_ ) ).at(m,n), ( reflhs_ * eval( refrhs_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * orhs_ )(m,n), ( eval( reflhs_ ) * refrhs_ )(m,n) ) ||
          !equal( ( eval( lhs_ ) * orhs_ ).at(m,n), ( eval( reflhs_ ) * refrhs_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( lhs_ ) * eval( orhs_ ) )(m,n), ( eval( reflhs_ ) * eval( refrhs_ ) )(m,n) ) ||
          !equal( ( eval( lhs_ ) * eval( orhs_ ) ).at(m,n), ( eval( reflhs_ ) * eval( refrhs_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT1 ).name() << "\n"
             << "   Right-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      ( lhs_ * orhs_ ).at( 0UL, orhs_.columns() );

      std::ostringstream oss;
      oss << " Test : Checked element access of multiplication expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side row-major sparse matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Right-hand side column-major sparse matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      ( lhs_ * orhs_ ).at( lhs_.rows(), 0UL );

      std::ostringstream oss;
      oss << " Test : Checked element access of multiplication expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side row-major sparse matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   Right-hand side column-major sparse matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}


   //=====================================================================================
   // Testing the element access with a column-major matrix and a row-major matrix
   //=====================================================================================

   if( olhs_.rows() > 0UL && rhs_.columns() > 0UL )
   {
      const size_t m( olhs_.rows()   - 1UL );
      const size_t n( rhs_.columns() - 1UL );

      if( !equal( ( olhs_ * rhs_ )(m,n), ( reflhs_ * refrhs_ )(m,n) ) ||
          !equal( ( olhs_ * rhs_ ).at(m,n), ( reflhs_ * refrhs_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( olhs_ * eval( rhs_ ) )(m,n), ( reflhs_ * eval( refrhs_ ) )(m,n) ) ||
          !equal( ( olhs_ * eval( rhs_ ) ).at(m,n), ( reflhs_ * eval( refrhs_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( olhs_ ) * rhs_ )(m,n), ( eval( reflhs_ ) * refrhs_ )(m,n) ) ||
          !equal( ( eval( olhs_ ) * rhs_ ).at(m,n), ( eval( reflhs_ ) * refrhs_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( olhs_ ) * eval( rhs_ ) )(m,n), ( eval( reflhs_ ) * eval( refrhs_ ) )(m,n) ) ||
          !equal( ( eval( olhs_ ) * eval( rhs_ ) ).at(m,n), ( eval( reflhs_ ) * eval( refrhs_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side row-major sparse matrix type:\n"
             << "     " << typeid( MT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      ( olhs_ * rhs_ ).at( 0UL, rhs_.columns() );

      std::ostringstream oss;
      oss << " Test : Checked element access of multiplication expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side column-major sparse matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Right-hand side row-major sparse matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      ( olhs_ * rhs_ ).at( olhs_.rows(), 0UL );

      std::ostringstream oss;
      oss << " Test : Checked element access of multiplication expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side column-major sparse matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Right-hand side row-major sparse matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}


   //=====================================================================================
   // Testing the element access with two column-major matrices
   //=====================================================================================

   if( olhs_.rows() > 0UL && orhs_.columns() > 0UL )
   {
      const size_t m( olhs_.rows()    - 1UL );
      const size_t n( orhs_.columns() - 1UL );

      if( !equal( ( olhs_ * orhs_ )(m,n), ( reflhs_ * refrhs_ )(m,n) ) ||
          !equal( ( olhs_ * orhs_ ).at(m,n), ( reflhs_ * refrhs_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( olhs_ * eval( orhs_ ) )(m,n), ( reflhs_ * eval( refrhs_ ) )(m,n) ) ||
          !equal( ( olhs_ * eval( orhs_ ) ).at(m,n), ( reflhs_ * eval( refrhs_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( olhs_ ) * orhs_ )(m,n), ( eval( reflhs_ ) * refrhs_ )(m,n) ) ||
          !equal( ( eval( olhs_ ) * orhs_ ).at(m,n), ( eval( reflhs_ ) * refrhs_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( ( eval( olhs_ ) * eval( orhs_ ) )(m,n), ( eval( reflhs_ ) * eval( refrhs_ ) )(m,n) ) ||
          !equal( ( eval( olhs_ ) * eval( orhs_ ) ).at(m,n), ( eval( reflhs_ ) * eval( refrhs_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated multiplication expression\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT1 ).name() << "\n"
             << "   Right-hand side column-major sparse matrix type:\n"
             << "     " << typeid( OMT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      ( olhs_ * orhs_ ).at( 0UL, orhs_.columns() );

      std::ostringstream oss;
      oss << " Test : Checked element access of multiplication expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side column-major sparse matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Right-hand side column-major sparse matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      ( olhs_ * orhs_ ).at( olhs_.rows(), 0UL );

      std::ostringstream oss;
      oss << " Test : Checked element access of multiplication expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side column-major sparse matrix type:\n"
          << "     " << typeid( OMT1 ).name() << "\n"
          << "   Right-hand side column-major sparse matrix type:\n"
          << "     " << typeid( OMT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the plain matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Multiplication
      //=====================================================================================

      // Multiplication with the given matrices
      {
         test_  = "Multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = lhs_ * rhs_;
            odres_  = lhs_ * rhs_;
            sres_   = lhs_ * rhs_;
            osres_  = lhs_ * rhs_;
            refres_ = reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = lhs_ * orhs_;
            odres_  = lhs_ * orhs_;
            sres_   = lhs_ * orhs_;
            osres_  = lhs_ * orhs_;
            refres_ = reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = olhs_ * rhs_;
            odres_  = olhs_ * rhs_;
            sres_   = olhs_ * rhs_;
            osres_  = olhs_ * rhs_;
            refres_ = reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = olhs_ * orhs_;
            odres_  = olhs_ * orhs_;
            sres_   = olhs_ * orhs_;
            osres_  = olhs_ * orhs_;
            refres_ = reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Multiplication with evaluated matrices
      {
         test_  = "Multiplication with evaluated matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = eval( lhs_ ) * eval( rhs_ );
            odres_  = eval( lhs_ ) * eval( rhs_ );
            sres_   = eval( lhs_ ) * eval( rhs_ );
            osres_  = eval( lhs_ ) * eval( rhs_ );
            refres_ = eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = eval( lhs_ ) * eval( orhs_ );
            odres_  = eval( lhs_ ) * eval( orhs_ );
            sres_   = eval( lhs_ ) * eval( orhs_ );
            osres_  = eval( lhs_ ) * eval( orhs_ );
            refres_ = eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = eval( olhs_ ) * eval( rhs_ );
            odres_  = eval( olhs_ ) * eval( rhs_ );
            sres_   = eval( olhs_ ) * eval( rhs_ );
            osres_  = eval( olhs_ ) * eval( rhs_ );
            refres_ = eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = eval( olhs_ ) * eval( orhs_ );
            odres_  = eval( olhs_ ) * eval( orhs_ );
            sres_   = eval( olhs_ ) * eval( orhs_ );
            osres_  = eval( olhs_ ) * eval( orhs_ );
            refres_ = eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Multiplication with addition assignment
      //=====================================================================================

      // Multiplication with addition assignment with the given matrices
      {
         test_  = "Multiplication with addition assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += lhs_ * orhs_;
            odres_  += lhs_ * orhs_;
            sres_   += lhs_ * orhs_;
            osres_  += lhs_ * orhs_;
            refres_ += reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += olhs_ * rhs_;
            odres_  += olhs_ * rhs_;
            sres_   += olhs_ * rhs_;
            osres_  += olhs_ * rhs_;
            refres_ += reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += olhs_ * orhs_;
            odres_  += olhs_ * orhs_;
            sres_   += olhs_ * orhs_;
            osres_  += olhs_ * orhs_;
            refres_ += reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Multiplication with addition assignment with evaluated matrices
      {
         test_  = "Multiplication with addition assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += eval( lhs_ ) * eval( orhs_ );
            odres_  += eval( lhs_ ) * eval( orhs_ );
            sres_   += eval( lhs_ ) * eval( orhs_ );
            osres_  += eval( lhs_ ) * eval( orhs_ );
            refres_ += eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += eval( olhs_ ) * eval( rhs_ );
            odres_  += eval( olhs_ ) * eval( rhs_ );
            sres_   += eval( olhs_ ) * eval( rhs_ );
            osres_  += eval( olhs_ ) * eval( rhs_ );
            refres_ += eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += eval( olhs_ ) * eval( orhs_ );
            odres_  += eval( olhs_ ) * eval( orhs_ );
            sres_   += eval( olhs_ ) * eval( orhs_ );
            osres_  += eval( olhs_ ) * eval( orhs_ );
            refres_ += eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Multiplication with subtraction assignment with the given matrices
      //=====================================================================================

      // Multiplication with subtraction assignment with the given matrices
      {
         test_  = "Multiplication with subtraction assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= lhs_ * orhs_;
            odres_  -= lhs_ * orhs_;
            sres_   -= lhs_ * orhs_;
            osres_  -= lhs_ * orhs_;
            refres_ -= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= olhs_ * rhs_;
            odres_  -= olhs_ * rhs_;
            sres_   -= olhs_ * rhs_;
            osres_  -= olhs_ * rhs_;
            refres_ -= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= olhs_ * orhs_;
            odres_  -= olhs_ * orhs_;
            sres_   -= olhs_ * orhs_;
            osres_  -= olhs_ * orhs_;
            refres_ -= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Multiplication with subtraction assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= eval( lhs_ ) * eval( orhs_ );
            odres_  -= eval( lhs_ ) * eval( orhs_ );
            sres_   -= eval( lhs_ ) * eval( orhs_ );
            osres_  -= eval( lhs_ ) * eval( orhs_ );
            refres_ -= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= eval( olhs_ ) * eval( rhs_ );
            odres_  -= eval( olhs_ ) * eval( rhs_ );
            sres_   -= eval( olhs_ ) * eval( rhs_ );
            osres_  -= eval( olhs_ ) * eval( rhs_ );
            refres_ -= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= eval( olhs_ ) * eval( orhs_ );
            odres_  -= eval( olhs_ ) * eval( orhs_ );
            sres_   -= eval( olhs_ ) * eval( orhs_ );
            osres_  -= eval( olhs_ ) * eval( orhs_ );
            refres_ -= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Multiplication with Schur product assignment
      //=====================================================================================

      // Multiplication with Schur product assignment with the given matrices
      {
         test_  = "Multiplication with Schur product assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= lhs_ * orhs_;
            odres_  %= lhs_ * orhs_;
            sres_   %= lhs_ * orhs_;
            osres_  %= lhs_ * orhs_;
            refres_ %= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= olhs_ * rhs_;
            odres_  %= olhs_ * rhs_;
            sres_   %= olhs_ * rhs_;
            osres_  %= olhs_ * rhs_;
            refres_ %= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= olhs_ * orhs_;
            odres_  %= olhs_ * orhs_;
            sres_   %= olhs_ * orhs_;
            osres_  %= olhs_ * orhs_;
            refres_ %= reflhs_ * refrhs_;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Multiplication with Schur product assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= eval( lhs_ ) * eval( orhs_ );
            odres_  %= eval( lhs_ ) * eval( orhs_ );
            sres_   %= eval( lhs_ ) * eval( orhs_ );
            osres_  %= eval( lhs_ ) * eval( orhs_ );
            refres_ %= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= eval( olhs_ ) * eval( rhs_ );
            odres_  %= eval( olhs_ ) * eval( rhs_ );
            sres_   %= eval( olhs_ ) * eval( rhs_ );
            osres_  %= eval( olhs_ ) * eval( rhs_ );
            refres_ %= eval( reflhs_ ) * eval( refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= eval( olhs_ ) * eval( orhs_ );
            odres_  %= eval( olhs_ ) * eval( orhs_ );
            sres_   %= eval( olhs_ ) * eval( orhs_ );
            osres_  %= eval( olhs_ ) * eval( orhs_ );
            refres_ %= eval( reflhs_ ) * eval( refrhs_ );
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
/*!\brief Testing the negated sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the negated matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      //=====================================================================================
      // Negated multiplication
      //=====================================================================================

      // Negated multiplication with the given matrices
      {
         test_  = "Negated multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = -( lhs_ * rhs_ );
            odres_  = -( lhs_ * rhs_ );
            sres_   = -( lhs_ * rhs_ );
            osres_  = -( lhs_ * rhs_ );
            refres_ = -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = -( lhs_ * orhs_ );
            odres_  = -( lhs_ * orhs_ );
            sres_   = -( lhs_ * orhs_ );
            osres_  = -( lhs_ * orhs_ );
            refres_ = -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = -( olhs_ * rhs_ );
            odres_  = -( olhs_ * rhs_ );
            sres_   = -( olhs_ * rhs_ );
            osres_  = -( olhs_ * rhs_ );
            refres_ = -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = -( olhs_ * orhs_ );
            odres_  = -( olhs_ * orhs_ );
            sres_   = -( olhs_ * orhs_ );
            osres_  = -( olhs_ * orhs_ );
            refres_ = -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Negated multiplication with evaluated matrices
      {
         test_  = "Negated multiplication with evaluated matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = -( eval( lhs_ ) * eval( rhs_ ) );
            odres_  = -( eval( lhs_ ) * eval( rhs_ ) );
            sres_   = -( eval( lhs_ ) * eval( rhs_ ) );
            osres_  = -( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = -( eval( lhs_ ) * eval( orhs_ ) );
            odres_  = -( eval( lhs_ ) * eval( orhs_ ) );
            sres_   = -( eval( lhs_ ) * eval( orhs_ ) );
            osres_  = -( eval( lhs_ ) * eval( orhs_ ) );
            refres_ = -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = -( eval( olhs_ ) * eval( rhs_ ) );
            odres_  = -( eval( olhs_ ) * eval( rhs_ ) );
            sres_   = -( eval( olhs_ ) * eval( rhs_ ) );
            osres_  = -( eval( olhs_ ) * eval( rhs_ ) );
            refres_ = -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = -( eval( olhs_ ) * eval( orhs_ ) );
            odres_  = -( eval( olhs_ ) * eval( orhs_ ) );
            sres_   = -( eval( olhs_ ) * eval( orhs_ ) );
            osres_  = -( eval( olhs_ ) * eval( orhs_ ) );
            refres_ = -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Negated multiplication with addition assignment
      //=====================================================================================

      // Negated multiplication with addition assignment with the given matrices
      {
         test_  = "Negated multiplication with addition assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += -( lhs_ * orhs_ );
            odres_  += -( lhs_ * orhs_ );
            sres_   += -( lhs_ * orhs_ );
            osres_  += -( lhs_ * orhs_ );
            refres_ += -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += -( olhs_ * rhs_ );
            odres_  += -( olhs_ * rhs_ );
            sres_   += -( olhs_ * rhs_ );
            osres_  += -( olhs_ * rhs_ );
            refres_ += -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += -( olhs_ * orhs_ );
            odres_  += -( olhs_ * orhs_ );
            sres_   += -( olhs_ * orhs_ );
            osres_  += -( olhs_ * orhs_ );
            refres_ += -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Negated multiplication with addition assignment with the given matrices
      {
         test_  = "Negated multiplication with addition assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += -( eval( lhs_ ) * eval( orhs_ ) );
            odres_  += -( eval( lhs_ ) * eval( orhs_ ) );
            sres_   += -( eval( lhs_ ) * eval( orhs_ ) );
            osres_  += -( eval( lhs_ ) * eval( orhs_ ) );
            refres_ += -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += -( eval( olhs_ ) * eval( rhs_ ) );
            odres_  += -( eval( olhs_ ) * eval( rhs_ ) );
            sres_   += -( eval( olhs_ ) * eval( rhs_ ) );
            osres_  += -( eval( olhs_ ) * eval( rhs_ ) );
            refres_ += -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += -( eval( olhs_ ) * eval( orhs_ ) );
            odres_  += -( eval( olhs_ ) * eval( orhs_ ) );
            sres_   += -( eval( olhs_ ) * eval( orhs_ ) );
            osres_  += -( eval( olhs_ ) * eval( orhs_ ) );
            refres_ += -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Negated multiplication with subtraction assignment
      //=====================================================================================

      // Negated multiplication with subtraction assignment with the given matrices
      {
         test_  = "Negated multiplication with subtraction assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= -( lhs_ * orhs_ );
            odres_  -= -( lhs_ * orhs_ );
            sres_   -= -( lhs_ * orhs_ );
            osres_  -= -( lhs_ * orhs_ );
            refres_ -= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= -( olhs_ * rhs_ );
            odres_  -= -( olhs_ * rhs_ );
            sres_   -= -( olhs_ * rhs_ );
            osres_  -= -( olhs_ * rhs_ );
            refres_ -= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= -( olhs_ * orhs_ );
            odres_  -= -( olhs_ * orhs_ );
            sres_   -= -( olhs_ * orhs_ );
            osres_  -= -( olhs_ * orhs_ );
            refres_ -= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Negated multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Negated multiplication with subtraction assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= -( eval( lhs_ ) * eval( orhs_ ) );
            odres_  -= -( eval( lhs_ ) * eval( orhs_ ) );
            sres_   -= -( eval( lhs_ ) * eval( orhs_ ) );
            osres_  -= -( eval( lhs_ ) * eval( orhs_ ) );
            refres_ -= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= -( eval( olhs_ ) * eval( rhs_ ) );
            odres_  -= -( eval( olhs_ ) * eval( rhs_ ) );
            sres_   -= -( eval( olhs_ ) * eval( rhs_ ) );
            osres_  -= -( eval( olhs_ ) * eval( rhs_ ) );
            refres_ -= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= -( eval( olhs_ ) * eval( orhs_ ) );
            odres_  -= -( eval( olhs_ ) * eval( orhs_ ) );
            sres_   -= -( eval( olhs_ ) * eval( orhs_ ) );
            osres_  -= -( eval( olhs_ ) * eval( orhs_ ) );
            refres_ -= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Negated multiplication with Schur product assignment
      //=====================================================================================

      // Negated multiplication with Schur product assignment with the given matrices
      {
         test_  = "Negated multiplication with Schur product assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= -( lhs_ * orhs_ );
            odres_  %= -( lhs_ * orhs_ );
            sres_   %= -( lhs_ * orhs_ );
            osres_  %= -( lhs_ * orhs_ );
            refres_ %= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= -( olhs_ * rhs_ );
            odres_  %= -( olhs_ * rhs_ );
            sres_   %= -( olhs_ * rhs_ );
            osres_  %= -( olhs_ * rhs_ );
            refres_ %= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= -( olhs_ * orhs_ );
            odres_  %= -( olhs_ * orhs_ );
            sres_   %= -( olhs_ * orhs_ );
            osres_  %= -( olhs_ * orhs_ );
            refres_ %= -( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Negated multiplication with Schur product assignment with the given matrices
      {
         test_  = "Negated multiplication with Schur product assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= -( eval( lhs_ ) * eval( orhs_ ) );
            odres_  %= -( eval( lhs_ ) * eval( orhs_ ) );
            sres_   %= -( eval( lhs_ ) * eval( orhs_ ) );
            osres_  %= -( eval( lhs_ ) * eval( orhs_ ) );
            refres_ %= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= -( eval( olhs_ ) * eval( rhs_ ) );
            odres_  %= -( eval( olhs_ ) * eval( rhs_ ) );
            sres_   %= -( eval( olhs_ ) * eval( rhs_ ) );
            osres_  %= -( eval( olhs_ ) * eval( rhs_ ) );
            refres_ %= -( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= -( eval( olhs_ ) * eval( orhs_ ) );
            odres_  %= -( eval( olhs_ ) * eval( orhs_ ) );
            sres_   %= -( eval( olhs_ ) * eval( orhs_ ) );
            osres_  %= -( eval( olhs_ ) * eval( orhs_ ) );
            refres_ %= -( eval( reflhs_ ) * eval( refrhs_ ) );
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
/*!\brief Testing the scaled sparse matrix/sparse matrix multiplication.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the scaled matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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

         checkResults<MT1,MT2>();
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

         checkResults<MT1,MT2>();
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

         checkResults<MT1,MT2>();
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

         checkResults<MT1,MT2>();
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

         checkResults<MT1,MT2>();
      }


      //=====================================================================================
      // Scaled multiplication (s*OP)
      //=====================================================================================

      // Scaled multiplication with the given matrices
      {
         test_  = "Scaled multiplication with the given matrices (s*OP)";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = scalar * ( lhs_ * rhs_ );
            odres_  = scalar * ( lhs_ * rhs_ );
            sres_   = scalar * ( lhs_ * rhs_ );
            osres_  = scalar * ( lhs_ * rhs_ );
            refres_ = scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = scalar * ( lhs_ * orhs_ );
            odres_  = scalar * ( lhs_ * orhs_ );
            sres_   = scalar * ( lhs_ * orhs_ );
            osres_  = scalar * ( lhs_ * orhs_ );
            refres_ = scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = scalar * ( olhs_ * rhs_ );
            odres_  = scalar * ( olhs_ * rhs_ );
            sres_   = scalar * ( olhs_ * rhs_ );
            osres_  = scalar * ( olhs_ * rhs_ );
            refres_ = scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = scalar * ( olhs_ * orhs_ );
            odres_  = scalar * ( olhs_ * orhs_ );
            sres_   = scalar * ( olhs_ * orhs_ );
            osres_  = scalar * ( olhs_ * orhs_ );
            refres_ = scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled multiplication with evaluated matrices
      {
         test_  = "Scaled multiplication with evaluated matrices (s*OP)";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            odres_  = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            sres_   = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            osres_  = scalar * ( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            odres_  = scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            sres_   = scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            osres_  = scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            refres_ = scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            odres_  = scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            sres_   = scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            osres_  = scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            refres_ = scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            odres_  = scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            sres_   = scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            osres_  = scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            refres_ = scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled multiplication (OP*s)
      //=====================================================================================

      // Scaled multiplication with the given matrices
      {
         test_  = "Scaled multiplication with the given matrices (OP*s)";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = ( lhs_ * rhs_ ) * scalar;
            odres_  = ( lhs_ * rhs_ ) * scalar;
            sres_   = ( lhs_ * rhs_ ) * scalar;
            osres_  = ( lhs_ * rhs_ ) * scalar;
            refres_ = ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = ( lhs_ * orhs_ ) * scalar;
            odres_  = ( lhs_ * orhs_ ) * scalar;
            sres_   = ( lhs_ * orhs_ ) * scalar;
            osres_  = ( lhs_ * orhs_ ) * scalar;
            refres_ = ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = ( olhs_ * rhs_ ) * scalar;
            odres_  = ( olhs_ * rhs_ ) * scalar;
            sres_   = ( olhs_ * rhs_ ) * scalar;
            osres_  = ( olhs_ * rhs_ ) * scalar;
            refres_ = ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = ( olhs_ * orhs_ ) * scalar;
            odres_  = ( olhs_ * orhs_ ) * scalar;
            sres_   = ( olhs_ * orhs_ ) * scalar;
            osres_  = ( olhs_ * orhs_ ) * scalar;
            refres_ = ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled multiplication with evaluated matrices
      {
         test_  = "Scaled multiplication with evaluated matrices (OP*s)";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            odres_  = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            osres_  = ( eval( lhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            odres_  = ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            sres_   = ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            osres_  = ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            odres_  = ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   = ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            osres_  = ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            odres_  = ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            sres_   = ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            osres_  = ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled multiplication (OP/s)
      //=====================================================================================

      // Scaled multiplication with the given matrices
      {
         test_  = "Scaled multiplication with the given matrices (OP/s)";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = ( lhs_ * rhs_ ) / scalar;
            odres_  = ( lhs_ * rhs_ ) / scalar;
            sres_   = ( lhs_ * rhs_ ) / scalar;
            osres_  = ( lhs_ * rhs_ ) / scalar;
            refres_ = ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = ( lhs_ * orhs_ ) / scalar;
            odres_  = ( lhs_ * orhs_ ) / scalar;
            sres_   = ( lhs_ * orhs_ ) / scalar;
            osres_  = ( lhs_ * orhs_ ) / scalar;
            refres_ = ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = ( olhs_ * rhs_ ) / scalar;
            odres_  = ( olhs_ * rhs_ ) / scalar;
            sres_   = ( olhs_ * rhs_ ) / scalar;
            osres_  = ( olhs_ * rhs_ ) / scalar;
            refres_ = ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = ( olhs_ * orhs_ ) / scalar;
            odres_  = ( olhs_ * orhs_ ) / scalar;
            sres_   = ( olhs_ * orhs_ ) / scalar;
            osres_  = ( olhs_ * orhs_ ) / scalar;
            refres_ = ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled multiplication with evaluated matrices
      {
         test_  = "Scaled multiplication with evaluated matrices (OP/s)";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            odres_  = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            osres_  = ( eval( lhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            odres_  = ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            sres_   = ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            osres_  = ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            odres_  = ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   = ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            osres_  = ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            odres_  = ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            sres_   = ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            osres_  = ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            refres_ = ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled multiplication with addition assignment (s*OP)
      //=====================================================================================

      // Scaled multiplication with addition assignment with the given matrices
      {
         test_  = "Scaled multiplication with addition assignment with the given matrices (s*OP)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += scalar * ( lhs_ * orhs_ );
            odres_  += scalar * ( lhs_ * orhs_ );
            sres_   += scalar * ( lhs_ * orhs_ );
            osres_  += scalar * ( lhs_ * orhs_ );
            refres_ += scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += scalar * ( olhs_ * rhs_ );
            odres_  += scalar * ( olhs_ * rhs_ );
            sres_   += scalar * ( olhs_ * rhs_ );
            osres_  += scalar * ( olhs_ * rhs_ );
            refres_ += scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += scalar * ( olhs_ * orhs_ );
            odres_  += scalar * ( olhs_ * orhs_ );
            sres_   += scalar * ( olhs_ * orhs_ );
            osres_  += scalar * ( olhs_ * orhs_ );
            refres_ += scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled multiplication with addition assignment with evaluated matrices
      {
         test_  = "Scaled multiplication with addition assignment with evaluated matrices (s*OP)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            odres_  += scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            sres_   += scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            osres_  += scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            odres_  += scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            sres_   += scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            osres_  += scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            odres_  += scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            sres_   += scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            osres_  += scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            refres_ += scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled multiplication with addition assignment (OP*s)
      //=====================================================================================

      // Scaled multiplication with addition assignment with the given matrices
      {
         test_  = "Scaled multiplication with addition assignment with the given matrices (OP*s)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += ( lhs_ * orhs_ ) * scalar;
            odres_  += ( lhs_ * orhs_ ) * scalar;
            sres_   += ( lhs_ * orhs_ ) * scalar;
            osres_  += ( lhs_ * orhs_ ) * scalar;
            refres_ += ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += ( olhs_ * rhs_ ) * scalar;
            odres_  += ( olhs_ * rhs_ ) * scalar;
            sres_   += ( olhs_ * rhs_ ) * scalar;
            osres_  += ( olhs_ * rhs_ ) * scalar;
            refres_ += ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += ( olhs_ * orhs_ ) * scalar;
            odres_  += ( olhs_ * orhs_ ) * scalar;
            sres_   += ( olhs_ * orhs_ ) * scalar;
            osres_  += ( olhs_ * orhs_ ) * scalar;
            refres_ += ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled multiplication with addition assignment with evaluated matrices
      {
         test_  = "Scaled multiplication with addition assignment with evaluated matrices (OP*s)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            odres_  += ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            sres_   += ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            osres_  += ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            odres_  += ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   += ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            osres_  += ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            odres_  += ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            sres_   += ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            osres_  += ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled multiplication with addition assignment (OP/s)
      //=====================================================================================

      // Scaled multiplication with addition assignment with the given matrices
      {
         test_  = "Scaled multiplication with addition assignment with the given matrices (OP/s)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += ( lhs_ * orhs_ ) / scalar;
            odres_  += ( lhs_ * orhs_ ) / scalar;
            sres_   += ( lhs_ * orhs_ ) / scalar;
            osres_  += ( lhs_ * orhs_ ) / scalar;
            refres_ += ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += ( olhs_ * rhs_ ) / scalar;
            odres_  += ( olhs_ * rhs_ ) / scalar;
            sres_   += ( olhs_ * rhs_ ) / scalar;
            osres_  += ( olhs_ * rhs_ ) / scalar;
            refres_ += ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += ( olhs_ * orhs_ ) / scalar;
            odres_  += ( olhs_ * orhs_ ) / scalar;
            sres_   += ( olhs_ * orhs_ ) / scalar;
            osres_  += ( olhs_ * orhs_ ) / scalar;
            refres_ += ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled multiplication with addition assignment with evaluated matrices
      {
         test_  = "Scaled multiplication with addition assignment with evaluated matrices (OP/s)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            odres_  += ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            sres_   += ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            osres_  += ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            odres_  += ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   += ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            osres_  += ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            odres_  += ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            sres_   += ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            osres_  += ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            refres_ += ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled multiplication with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled multiplication with subtraction assignment with the given matrices
      {
         test_  = "Scaled multiplication with subtraction assignment with the given matrices (s*OP)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= scalar * ( lhs_ * orhs_ );
            odres_  -= scalar * ( lhs_ * orhs_ );
            sres_   -= scalar * ( lhs_ * orhs_ );
            osres_  -= scalar * ( lhs_ * orhs_ );
            refres_ -= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= scalar * ( olhs_ * rhs_ );
            odres_  -= scalar * ( olhs_ * rhs_ );
            sres_   -= scalar * ( olhs_ * rhs_ );
            osres_  -= scalar * ( olhs_ * rhs_ );
            refres_ -= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= scalar * ( olhs_ * orhs_ );
            odres_  -= scalar * ( olhs_ * orhs_ );
            sres_   -= scalar * ( olhs_ * orhs_ );
            osres_  -= scalar * ( olhs_ * orhs_ );
            refres_ -= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Scaled multiplication with subtraction assignment with evaluated matrices (s*OP)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            odres_  -= scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            sres_   -= scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            osres_  -= scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            odres_  -= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            sres_   -= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            osres_  -= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            odres_  -= scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            sres_   -= scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            osres_  -= scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            refres_ -= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled multiplication with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled multiplication with subtraction assignment with the given matrices
      {
         test_  = "Scaled multiplication with subtraction assignment with the given matrices (OP*s)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= ( lhs_ * orhs_ ) * scalar;
            odres_  -= ( lhs_ * orhs_ ) * scalar;
            sres_   -= ( lhs_ * orhs_ ) * scalar;
            osres_  -= ( lhs_ * orhs_ ) * scalar;
            refres_ -= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= ( olhs_ * rhs_ ) * scalar;
            odres_  -= ( olhs_ * rhs_ ) * scalar;
            sres_   -= ( olhs_ * rhs_ ) * scalar;
            osres_  -= ( olhs_ * rhs_ ) * scalar;
            refres_ -= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= ( olhs_ * orhs_ ) * scalar;
            odres_  -= ( olhs_ * orhs_ ) * scalar;
            sres_   -= ( olhs_ * orhs_ ) * scalar;
            osres_  -= ( olhs_ * orhs_ ) * scalar;
            refres_ -= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Scaled multiplication with subtraction assignment with evaluated matrices (OP*s)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            odres_  -= ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            sres_   -= ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            osres_  -= ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            odres_  -= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   -= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            osres_  -= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            odres_  -= ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            sres_   -= ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            osres_  -= ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled multiplication with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled multiplication with subtraction assignment with the given matrices
      {
         test_  = "Scaled multiplication with subtraction assignment with the given matrices (OP/s)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= ( lhs_ * orhs_ ) / scalar;
            odres_  -= ( lhs_ * orhs_ ) / scalar;
            sres_   -= ( lhs_ * orhs_ ) / scalar;
            osres_  -= ( lhs_ * orhs_ ) / scalar;
            refres_ -= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= ( olhs_ * rhs_ ) / scalar;
            odres_  -= ( olhs_ * rhs_ ) / scalar;
            sres_   -= ( olhs_ * rhs_ ) / scalar;
            osres_  -= ( olhs_ * rhs_ ) / scalar;
            refres_ -= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= ( olhs_ * orhs_ ) / scalar;
            odres_  -= ( olhs_ * orhs_ ) / scalar;
            sres_   -= ( olhs_ * orhs_ ) / scalar;
            osres_  -= ( olhs_ * orhs_ ) / scalar;
            refres_ -= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Scaled multiplication with subtraction assignment with evaluated matrices (OP/s)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            odres_  -= ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            sres_   -= ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            osres_  -= ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            odres_  -= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   -= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            osres_  -= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            odres_  -= ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            sres_   -= ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            osres_  -= ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            refres_ -= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled multiplication with Schur product assignment (s*OP)
      //=====================================================================================

      // Scaled multiplication with Schur product assignment with the given matrices
      {
         test_  = "Scaled multiplication with Schur product assignment with the given matrices (s*OP)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= scalar * ( lhs_ * orhs_ );
            odres_  %= scalar * ( lhs_ * orhs_ );
            sres_   %= scalar * ( lhs_ * orhs_ );
            osres_  %= scalar * ( lhs_ * orhs_ );
            refres_ %= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= scalar * ( olhs_ * rhs_ );
            odres_  %= scalar * ( olhs_ * rhs_ );
            sres_   %= scalar * ( olhs_ * rhs_ );
            osres_  %= scalar * ( olhs_ * rhs_ );
            refres_ %= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= scalar * ( olhs_ * orhs_ );
            odres_  %= scalar * ( olhs_ * orhs_ );
            sres_   %= scalar * ( olhs_ * orhs_ );
            osres_  %= scalar * ( olhs_ * orhs_ );
            refres_ %= scalar * ( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Scaled multiplication with Schur product assignment with evaluated matrices (s*OP)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            odres_  %= scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            sres_   %= scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            osres_  %= scalar * ( eval( lhs_ ) * eval( orhs_ ) );
            refres_ %= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            odres_  %= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            sres_   %= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            osres_  %= scalar * ( eval( olhs_ ) * eval( rhs_ ) );
            refres_ %= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            odres_  %= scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            sres_   %= scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            osres_  %= scalar * ( eval( olhs_ ) * eval( orhs_ ) );
            refres_ %= scalar * ( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled multiplication with Schur product assignment (OP*s)
      //=====================================================================================

      // Scaled multiplication with Schur product assignment with the given matrices
      {
         test_  = "Scaled multiplication with Schur product assignment with the given matrices (OP*s)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= ( lhs_ * orhs_ ) * scalar;
            odres_  %= ( lhs_ * orhs_ ) * scalar;
            sres_   %= ( lhs_ * orhs_ ) * scalar;
            osres_  %= ( lhs_ * orhs_ ) * scalar;
            refres_ %= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= ( olhs_ * rhs_ ) * scalar;
            odres_  %= ( olhs_ * rhs_ ) * scalar;
            sres_   %= ( olhs_ * rhs_ ) * scalar;
            osres_  %= ( olhs_ * rhs_ ) * scalar;
            refres_ %= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= ( olhs_ * orhs_ ) * scalar;
            odres_  %= ( olhs_ * orhs_ ) * scalar;
            sres_   %= ( olhs_ * orhs_ ) * scalar;
            osres_  %= ( olhs_ * orhs_ ) * scalar;
            refres_ %= ( reflhs_ * refrhs_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Scaled multiplication with Schur product assignment with evaluated matrices (OP*s)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            odres_  %= ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            sres_   %= ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            osres_  %= ( eval( lhs_ ) * eval( orhs_ ) ) * scalar;
            refres_ %= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            odres_  %= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            sres_   %= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            osres_  %= ( eval( olhs_ ) * eval( rhs_ ) ) * scalar;
            refres_ %= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            odres_  %= ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            sres_   %= ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            osres_  %= ( eval( olhs_ ) * eval( orhs_ ) ) * scalar;
            refres_ %= ( eval( reflhs_ ) * eval( refrhs_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Scaled multiplication with Schur product assignment (OP/s)
      //=====================================================================================

      // Scaled multiplication with Schur product assignment with the given matrices
      {
         test_  = "Scaled multiplication with Schur product assignment with the given matrices (OP/s)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= ( lhs_ * orhs_ ) / scalar;
            odres_  %= ( lhs_ * orhs_ ) / scalar;
            sres_   %= ( lhs_ * orhs_ ) / scalar;
            osres_  %= ( lhs_ * orhs_ ) / scalar;
            refres_ %= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= ( olhs_ * rhs_ ) / scalar;
            odres_  %= ( olhs_ * rhs_ ) / scalar;
            sres_   %= ( olhs_ * rhs_ ) / scalar;
            osres_  %= ( olhs_ * rhs_ ) / scalar;
            refres_ %= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= ( olhs_ * orhs_ ) / scalar;
            odres_  %= ( olhs_ * orhs_ ) / scalar;
            sres_   %= ( olhs_ * orhs_ ) / scalar;
            osres_  %= ( olhs_ * orhs_ ) / scalar;
            refres_ %= ( reflhs_ * refrhs_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Scaled multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Scaled multiplication with Schur product assignment with evaluated matrices (OP/s)";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            odres_  %= ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            sres_   %= ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            osres_  %= ( eval( lhs_ ) * eval( orhs_ ) ) / scalar;
            refres_ %= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            odres_  %= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            sres_   %= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            osres_  %= ( eval( olhs_ ) * eval( rhs_ ) ) / scalar;
            refres_ %= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            odres_  %= ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            sres_   %= ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            osres_  %= ( eval( olhs_ ) * eval( orhs_ ) ) / scalar;
            refres_ %= ( eval( reflhs_ ) * eval( refrhs_ ) ) / scalar;
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
/*!\brief Testing the transpose sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the transpose matrix multiplication with plain assignment. In case
// any error resulting from the multiplication or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANS_OPERATION > 1 )
   {
      //=====================================================================================
      // Transpose multiplication
      //=====================================================================================

      // Transpose multiplication with the given matrices
      {
         test_  = "Transpose multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initTransposeResults();
            tdres_  = trans( lhs_ * rhs_ );
            todres_ = trans( lhs_ * rhs_ );
            tsres_  = trans( lhs_ * rhs_ );
            tosres_ = trans( lhs_ * rhs_ );
            refres_ = trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkTransposeResults<MT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = trans( lhs_ * orhs_ );
            todres_ = trans( lhs_ * orhs_ );
            tsres_  = trans( lhs_ * orhs_ );
            tosres_ = trans( lhs_ * orhs_ );
            refres_ = trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkTransposeResults<MT1,OMT2>();

         try {
            initTransposeResults();
            tdres_  = trans( olhs_ * rhs_ );
            todres_ = trans( olhs_ * rhs_ );
            tsres_  = trans( olhs_ * rhs_ );
            tosres_ = trans( olhs_ * rhs_ );
            refres_ = trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkTransposeResults<OMT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = trans( olhs_ * orhs_ );
            todres_ = trans( olhs_ * orhs_ );
            tsres_  = trans( olhs_ * orhs_ );
            tosres_ = trans( olhs_ * orhs_ );
            refres_ = trans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkTransposeResults<OMT1,OMT2>();
      }

      // Transpose multiplication with evaluated matrices
      {
         test_  = "Transpose multiplication with evaluated matrices";
         error_ = "Failed multiplication operation";

         try {
            initTransposeResults();
            tdres_  = trans( eval( lhs_ ) * eval( rhs_ ) );
            todres_ = trans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_  = trans( eval( lhs_ ) * eval( rhs_ ) );
            tosres_ = trans( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkTransposeResults<MT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = trans( eval( lhs_ ) * eval( orhs_ ) );
            todres_ = trans( eval( lhs_ ) * eval( orhs_ ) );
            tsres_  = trans( eval( lhs_ ) * eval( orhs_ ) );
            tosres_ = trans( eval( lhs_ ) * eval( orhs_ ) );
            refres_ = trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkTransposeResults<MT1,OMT2>();

         try {
            initTransposeResults();
            tdres_  = trans( eval( olhs_ ) * eval( rhs_ ) );
            todres_ = trans( eval( olhs_ ) * eval( rhs_ ) );
            tsres_  = trans( eval( olhs_ ) * eval( rhs_ ) );
            tosres_ = trans( eval( olhs_ ) * eval( rhs_ ) );
            refres_ = trans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkTransposeResults<OMT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = trans( eval( olhs_ ) * eval( orhs_ ) );
            todres_ = trans( eval( olhs_ ) * eval( orhs_ ) );
            tsres_  = trans( eval( olhs_ ) * eval( orhs_ ) );
            tosres_ = trans( eval( olhs_ ) * eval( orhs_ ) );
            refres_ = trans( eval( reflhs_ ) * eval( refrhs_ ) );
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
/*!\brief Testing the conjugate transpose sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the conjugate transpose matrix multiplication with plain assignment. In
// case any error resulting from the multiplication or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testCTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION > 1 )
   {
      //=====================================================================================
      // Conjugate transpose multiplication
      //=====================================================================================

      // Conjugate transpose multiplication with the given matrices
      {
         test_  = "Conjugate transpose multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initTransposeResults();
            tdres_  = ctrans( lhs_ * rhs_ );
            todres_ = ctrans( lhs_ * rhs_ );
            tsres_  = ctrans( lhs_ * rhs_ );
            tosres_ = ctrans( lhs_ * rhs_ );
            refres_ = ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkTransposeResults<MT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = ctrans( lhs_ * orhs_ );
            todres_ = ctrans( lhs_ * orhs_ );
            tsres_  = ctrans( lhs_ * orhs_ );
            tosres_ = ctrans( lhs_ * orhs_ );
            refres_ = ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkTransposeResults<MT1,OMT2>();

         try {
            initTransposeResults();
            tdres_  = ctrans( olhs_ * rhs_ );
            todres_ = ctrans( olhs_ * rhs_ );
            tsres_  = ctrans( olhs_ * rhs_ );
            tosres_ = ctrans( olhs_ * rhs_ );
            refres_ = ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkTransposeResults<OMT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = ctrans( olhs_ * orhs_ );
            todres_ = ctrans( olhs_ * orhs_ );
            tsres_  = ctrans( olhs_ * orhs_ );
            tosres_ = ctrans( olhs_ * orhs_ );
            refres_ = ctrans( reflhs_ * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkTransposeResults<OMT1,OMT2>();
      }

      // Conjugate transpose multiplication with evaluated matrices
      {
         test_  = "Conjugate transpose multiplication with evaluated matrices";
         error_ = "Failed multiplication operation";

         try {
            initTransposeResults();
            tdres_  = ctrans( eval( lhs_ ) * eval( rhs_ ) );
            todres_ = ctrans( eval( lhs_ ) * eval( rhs_ ) );
            tsres_  = ctrans( eval( lhs_ ) * eval( rhs_ ) );
            tosres_ = ctrans( eval( lhs_ ) * eval( rhs_ ) );
            refres_ = ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkTransposeResults<MT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = ctrans( eval( lhs_ ) * eval( orhs_ ) );
            todres_ = ctrans( eval( lhs_ ) * eval( orhs_ ) );
            tsres_  = ctrans( eval( lhs_ ) * eval( orhs_ ) );
            tosres_ = ctrans( eval( lhs_ ) * eval( orhs_ ) );
            refres_ = ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkTransposeResults<MT1,OMT2>();

         try {
            initTransposeResults();
            tdres_  = ctrans( eval( olhs_ ) * eval( rhs_ ) );
            todres_ = ctrans( eval( olhs_ ) * eval( rhs_ ) );
            tsres_  = ctrans( eval( olhs_ ) * eval( rhs_ ) );
            tosres_ = ctrans( eval( olhs_ ) * eval( rhs_ ) );
            refres_ = ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkTransposeResults<OMT1,MT2>();

         try {
            initTransposeResults();
            tdres_  = ctrans( eval( olhs_ ) * eval( orhs_ ) );
            todres_ = ctrans( eval( olhs_ ) * eval( orhs_ ) );
            tsres_  = ctrans( eval( olhs_ ) * eval( orhs_ ) );
            tosres_ = ctrans( eval( olhs_ ) * eval( orhs_ ) );
            refres_ = ctrans( eval( reflhs_ ) * eval( refrhs_ ) );
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
/*!\brief Testing the abs sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the abs matrix multiplication with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// multiplication or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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
/*!\brief Testing the conjugate sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the conjugate matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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
/*!\brief Testing the \a real sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the \a real matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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
/*!\brief Testing the \a imag sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the \a imag matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testImagOperation()
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
/*!\brief Testing the evaluated sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the evaluated matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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
/*!\brief Testing the serialized sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the serialized matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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
/*!\brief Testing the non-aliased sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the non-aliased matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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
/*!\brief Testing the non-SIMD sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the non-SIMD matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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
/*!\brief Testing the symmetric sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the symmetric matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment in the case that the
// left-hand side element type is the lower data type. In case any error resulting from the
// multiplication or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testDeclSymOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION > 1 )
   {
      if( lhs_.rows() != rhs_.columns() )
         return;


      //=====================================================================================
      // Test-specific setup
      //=====================================================================================

      MT2  rhs;
      OMT2 orhs;
      RT2  refrhs;

      try {
         rhs    = trans( lhs_ );
         orhs   = trans( lhs_ );
         refrhs = trans( lhs_ );
      }
      catch( std::exception& ) {
         return;
      }


      //=====================================================================================
      // Declsym multiplication
      //=====================================================================================

      // Declsym multiplication with the given matrices
      {
         test_  = "Declsym multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = declsym( lhs_ * rhs );
            odres_  = declsym( lhs_ * rhs );
            sres_   = declsym( lhs_ * rhs );
            osres_  = declsym( lhs_ * rhs );
            refres_ = declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declsym( lhs_ * orhs );
            odres_  = declsym( lhs_ * orhs );
            sres_   = declsym( lhs_ * orhs );
            osres_  = declsym( lhs_ * orhs );
            refres_ = declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declsym( olhs_ * rhs );
            odres_  = declsym( olhs_ * rhs );
            sres_   = declsym( olhs_ * rhs );
            osres_  = declsym( olhs_ * rhs );
            refres_ = declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declsym( olhs_ * orhs );
            odres_  = declsym( olhs_ * orhs );
            sres_   = declsym( olhs_ * orhs );
            osres_  = declsym( olhs_ * orhs );
            refres_ = declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declsym multiplication with evaluated matrices
      {
         test_  = "Declsym multiplication with evaluated left-hand side matrix";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = declsym( eval( lhs_ ) * eval( rhs ) );
            odres_  = declsym( eval( lhs_ ) * eval( rhs ) );
            sres_   = declsym( eval( lhs_ ) * eval( rhs ) );
            osres_  = declsym( eval( lhs_ ) * eval( rhs ) );
            refres_ = declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declsym( eval( lhs_ ) * eval( orhs ) );
            odres_  = declsym( eval( lhs_ ) * eval( orhs ) );
            sres_   = declsym( eval( lhs_ ) * eval( orhs ) );
            osres_  = declsym( eval( lhs_ ) * eval( orhs ) );
            refres_ = declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declsym( eval( olhs_ ) * eval( rhs ) );
            odres_  = declsym( eval( olhs_ ) * eval( rhs ) );
            sres_   = declsym( eval( olhs_ ) * eval( rhs ) );
            osres_  = declsym( eval( olhs_ ) * eval( rhs ) );
            refres_ = declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declsym( eval( olhs_ ) * eval( orhs ) );
            odres_  = declsym( eval( olhs_ ) * eval( orhs ) );
            sres_   = declsym( eval( olhs_ ) * eval( orhs ) );
            osres_  = declsym( eval( olhs_ ) * eval( orhs ) );
            refres_ = declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declsym multiplication with addition assignment
      //=====================================================================================

      // Declsym multiplication with addition assignment with the given matrices
      {
         test_  = "Declsym multiplication with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declsym( lhs_ * rhs );
            odres_  += declsym( lhs_ * rhs );
            sres_   += declsym( lhs_ * rhs );
            osres_  += declsym( lhs_ * rhs );
            refres_ += declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declsym( lhs_ * orhs );
            odres_  += declsym( lhs_ * orhs );
            sres_   += declsym( lhs_ * orhs );
            osres_  += declsym( lhs_ * orhs );
            refres_ += declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declsym( olhs_ * rhs );
            odres_  += declsym( olhs_ * rhs );
            sres_   += declsym( olhs_ * rhs );
            osres_  += declsym( olhs_ * rhs );
            refres_ += declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declsym( olhs_ * orhs );
            odres_  += declsym( olhs_ * orhs );
            sres_   += declsym( olhs_ * orhs );
            osres_  += declsym( olhs_ * orhs );
            refres_ += declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declsym multiplication with addition assignment with evaluated matrices
      {
         test_  = "Declsym multiplication with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declsym( eval( lhs_ ) * eval( rhs ) );
            odres_  += declsym( eval( lhs_ ) * eval( rhs ) );
            sres_   += declsym( eval( lhs_ ) * eval( rhs ) );
            osres_  += declsym( eval( lhs_ ) * eval( rhs ) );
            refres_ += declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declsym( eval( lhs_ ) * eval( orhs ) );
            odres_  += declsym( eval( lhs_ ) * eval( orhs ) );
            sres_   += declsym( eval( lhs_ ) * eval( orhs ) );
            osres_  += declsym( eval( lhs_ ) * eval( orhs ) );
            refres_ += declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declsym( eval( olhs_ ) * eval( rhs ) );
            odres_  += declsym( eval( olhs_ ) * eval( rhs ) );
            sres_   += declsym( eval( olhs_ ) * eval( rhs ) );
            osres_  += declsym( eval( olhs_ ) * eval( rhs ) );
            refres_ += declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declsym( eval( olhs_ ) * eval( orhs ) );
            odres_  += declsym( eval( olhs_ ) * eval( orhs ) );
            sres_   += declsym( eval( olhs_ ) * eval( orhs ) );
            osres_  += declsym( eval( olhs_ ) * eval( orhs ) );
            refres_ += declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declsym multiplication with subtraction assignment
      //=====================================================================================

      // Declsym multiplication with subtraction assignment with the given matrices
      {
         test_  = "Declsym multiplication with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declsym( lhs_ * rhs );
            odres_  -= declsym( lhs_ * rhs );
            sres_   -= declsym( lhs_ * rhs );
            osres_  -= declsym( lhs_ * rhs );
            refres_ -= declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declsym( lhs_ * orhs );
            odres_  -= declsym( lhs_ * orhs );
            sres_   -= declsym( lhs_ * orhs );
            osres_  -= declsym( lhs_ * orhs );
            refres_ -= declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declsym( olhs_ * rhs );
            odres_  -= declsym( olhs_ * rhs );
            sres_   -= declsym( olhs_ * rhs );
            osres_  -= declsym( olhs_ * rhs );
            refres_ -= declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declsym( olhs_ * orhs );
            odres_  -= declsym( olhs_ * orhs );
            sres_   -= declsym( olhs_ * orhs );
            osres_  -= declsym( olhs_ * orhs );
            refres_ -= declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declsym multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Declsym multiplication with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declsym( eval( lhs_ ) * eval( rhs ) );
            odres_  -= declsym( eval( lhs_ ) * eval( rhs ) );
            sres_   -= declsym( eval( lhs_ ) * eval( rhs ) );
            osres_  -= declsym( eval( lhs_ ) * eval( rhs ) );
            refres_ -= declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declsym( eval( lhs_ ) * eval( orhs ) );
            odres_  -= declsym( eval( lhs_ ) * eval( orhs ) );
            sres_   -= declsym( eval( lhs_ ) * eval( orhs ) );
            osres_  -= declsym( eval( lhs_ ) * eval( orhs ) );
            refres_ -= declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declsym( eval( olhs_ ) * eval( rhs ) );
            odres_  -= declsym( eval( olhs_ ) * eval( rhs ) );
            sres_   -= declsym( eval( olhs_ ) * eval( rhs ) );
            osres_  -= declsym( eval( olhs_ ) * eval( rhs ) );
            refres_ -= declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declsym( eval( olhs_ ) * eval( orhs ) );
            odres_  -= declsym( eval( olhs_ ) * eval( orhs ) );
            sres_   -= declsym( eval( olhs_ ) * eval( orhs ) );
            osres_  -= declsym( eval( olhs_ ) * eval( orhs ) );
            refres_ -= declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declsym multiplication with Schur product assignment
      //=====================================================================================

      // Declsym multiplication with Schur product assignment with the given matrices
      {
         test_  = "Declsym multiplication with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declsym( lhs_ * rhs );
            odres_  %= declsym( lhs_ * rhs );
            sres_   %= declsym( lhs_ * rhs );
            osres_  %= declsym( lhs_ * rhs );
            refres_ %= declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declsym( lhs_ * orhs );
            odres_  %= declsym( lhs_ * orhs );
            sres_   %= declsym( lhs_ * orhs );
            osres_  %= declsym( lhs_ * orhs );
            refres_ %= declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declsym( olhs_ * rhs );
            odres_  %= declsym( olhs_ * rhs );
            sres_   %= declsym( olhs_ * rhs );
            osres_  %= declsym( olhs_ * rhs );
            refres_ %= declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declsym( olhs_ * orhs );
            odres_  %= declsym( olhs_ * orhs );
            sres_   %= declsym( olhs_ * orhs );
            osres_  %= declsym( olhs_ * orhs );
            refres_ %= declsym( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declsym multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Declsym multiplication with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declsym( eval( lhs_ ) * eval( rhs ) );
            odres_  %= declsym( eval( lhs_ ) * eval( rhs ) );
            sres_   %= declsym( eval( lhs_ ) * eval( rhs ) );
            osres_  %= declsym( eval( lhs_ ) * eval( rhs ) );
            refres_ %= declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declsym( eval( lhs_ ) * eval( orhs ) );
            odres_  %= declsym( eval( lhs_ ) * eval( orhs ) );
            sres_   %= declsym( eval( lhs_ ) * eval( orhs ) );
            osres_  %= declsym( eval( lhs_ ) * eval( orhs ) );
            refres_ %= declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declsym( eval( olhs_ ) * eval( rhs ) );
            odres_  %= declsym( eval( olhs_ ) * eval( rhs ) );
            sres_   %= declsym( eval( olhs_ ) * eval( rhs ) );
            osres_  %= declsym( eval( olhs_ ) * eval( rhs ) );
            refres_ %= declsym( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declsym( eval( olhs_ ) * eval( orhs ) );
            odres_  %= declsym( eval( olhs_ ) * eval( orhs ) );
            sres_   %= declsym( eval( olhs_ ) * eval( orhs ) );
            osres_  %= declsym( eval( olhs_ ) * eval( orhs ) );
            refres_ %= declsym( eval( reflhs_ ) * eval( refrhs ) );
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
/*!\brief Testing the symmetric sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the symmetric matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment in the case that the
// right-hand side element type is the lower data type. In case any error resulting from
// the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testDeclSymOperation( blaze::FalseType )
{
#if BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION > 1 )
   {
      if( lhs_.rows() != rhs_.columns() )
         return;


      //=====================================================================================
      // Test-specific setup
      //=====================================================================================

      MT1  lhs;
      OMT1 olhs;
      RT1  reflhs;

      try {
         lhs    = trans( rhs_ );
         olhs   = trans( rhs_ );
         reflhs = trans( rhs_ );
      }
      catch( std::exception& ) {
         return;
      }


      //=====================================================================================
      // Declsym multiplication
      //=====================================================================================

      // Declsym multiplication with the given matrices
      {
         test_  = "Declsym multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = declsym( lhs * rhs_ );
            odres_  = declsym( lhs * rhs_ );
            sres_   = declsym( lhs * rhs_ );
            osres_  = declsym( lhs * rhs_ );
            refres_ = declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declsym( lhs * orhs_ );
            odres_  = declsym( lhs * orhs_ );
            sres_   = declsym( lhs * orhs_ );
            osres_  = declsym( lhs * orhs_ );
            refres_ = declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declsym( olhs * rhs_ );
            odres_  = declsym( olhs * rhs_ );
            sres_   = declsym( olhs * rhs_ );
            osres_  = declsym( olhs * rhs_ );
            refres_ = declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declsym( olhs * orhs_ );
            odres_  = declsym( olhs * orhs_ );
            sres_   = declsym( olhs * orhs_ );
            osres_  = declsym( olhs * orhs_ );
            refres_ = declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declsym multiplication with evaluated matrices
      {
         test_  = "Declsym multiplication with evaluated left-hand side matrix";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = declsym( eval( lhs ) * eval( rhs_ ) );
            odres_  = declsym( eval( lhs ) * eval( rhs_ ) );
            sres_   = declsym( eval( lhs ) * eval( rhs_ ) );
            osres_  = declsym( eval( lhs ) * eval( rhs_ ) );
            refres_ = declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declsym( eval( lhs ) * eval( orhs_ ) );
            odres_  = declsym( eval( lhs ) * eval( orhs_ ) );
            sres_   = declsym( eval( lhs ) * eval( orhs_ ) );
            osres_  = declsym( eval( lhs ) * eval( orhs_ ) );
            refres_ = declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declsym( eval( olhs ) * eval( rhs_ ) );
            odres_  = declsym( eval( olhs ) * eval( rhs_ ) );
            sres_   = declsym( eval( olhs ) * eval( rhs_ ) );
            osres_  = declsym( eval( olhs ) * eval( rhs_ ) );
            refres_ = declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declsym( eval( olhs ) * eval( orhs_ ) );
            odres_  = declsym( eval( olhs ) * eval( orhs_ ) );
            sres_   = declsym( eval( olhs ) * eval( orhs_ ) );
            osres_  = declsym( eval( olhs ) * eval( orhs_ ) );
            refres_ = declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declsym multiplication with addition assignment
      //=====================================================================================

      // Declsym multiplication with addition assignment with the given matrices
      {
         test_  = "Declsym multiplication with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declsym( lhs * rhs_ );
            odres_  += declsym( lhs * rhs_ );
            sres_   += declsym( lhs * rhs_ );
            osres_  += declsym( lhs * rhs_ );
            refres_ += declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declsym( lhs * orhs_ );
            odres_  += declsym( lhs * orhs_ );
            sres_   += declsym( lhs * orhs_ );
            osres_  += declsym( lhs * orhs_ );
            refres_ += declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declsym( olhs * rhs_ );
            odres_  += declsym( olhs * rhs_ );
            sres_   += declsym( olhs * rhs_ );
            osres_  += declsym( olhs * rhs_ );
            refres_ += declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declsym( olhs * orhs_ );
            odres_  += declsym( olhs * orhs_ );
            sres_   += declsym( olhs * orhs_ );
            osres_  += declsym( olhs * orhs_ );
            refres_ += declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declsym multiplication with addition assignment with evaluated matrices
      {
         test_  = "Declsym multiplication with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declsym( eval( lhs ) * eval( rhs_ ) );
            odres_  += declsym( eval( lhs ) * eval( rhs_ ) );
            sres_   += declsym( eval( lhs ) * eval( rhs_ ) );
            osres_  += declsym( eval( lhs ) * eval( rhs_ ) );
            refres_ += declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declsym( eval( lhs ) * eval( orhs_ ) );
            odres_  += declsym( eval( lhs ) * eval( orhs_ ) );
            sres_   += declsym( eval( lhs ) * eval( orhs_ ) );
            osres_  += declsym( eval( lhs ) * eval( orhs_ ) );
            refres_ += declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declsym( eval( olhs ) * eval( rhs_ ) );
            odres_  += declsym( eval( olhs ) * eval( rhs_ ) );
            sres_   += declsym( eval( olhs ) * eval( rhs_ ) );
            osres_  += declsym( eval( olhs ) * eval( rhs_ ) );
            refres_ += declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declsym( eval( olhs ) * eval( orhs_ ) );
            odres_  += declsym( eval( olhs ) * eval( orhs_ ) );
            sres_   += declsym( eval( olhs ) * eval( orhs_ ) );
            osres_  += declsym( eval( olhs ) * eval( orhs_ ) );
            refres_ += declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declsym multiplication with subtraction assignment
      //=====================================================================================

      // Declsym multiplication with subtraction assignment with the given matrices
      {
         test_  = "Declsym multiplication with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declsym( lhs * rhs_ );
            odres_  -= declsym( lhs * rhs_ );
            sres_   -= declsym( lhs * rhs_ );
            osres_  -= declsym( lhs * rhs_ );
            refres_ -= declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declsym( lhs * orhs_ );
            odres_  -= declsym( lhs * orhs_ );
            sres_   -= declsym( lhs * orhs_ );
            osres_  -= declsym( lhs * orhs_ );
            refres_ -= declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declsym( olhs * rhs_ );
            odres_  -= declsym( olhs * rhs_ );
            sres_   -= declsym( olhs * rhs_ );
            osres_  -= declsym( olhs * rhs_ );
            refres_ -= declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declsym( olhs * orhs_ );
            odres_  -= declsym( olhs * orhs_ );
            sres_   -= declsym( olhs * orhs_ );
            osres_  -= declsym( olhs * orhs_ );
            refres_ -= declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declsym multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Declsym multiplication with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declsym( eval( lhs ) * eval( rhs_ ) );
            odres_  -= declsym( eval( lhs ) * eval( rhs_ ) );
            sres_   -= declsym( eval( lhs ) * eval( rhs_ ) );
            osres_  -= declsym( eval( lhs ) * eval( rhs_ ) );
            refres_ -= declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declsym( eval( lhs ) * eval( orhs_ ) );
            odres_  -= declsym( eval( lhs ) * eval( orhs_ ) );
            sres_   -= declsym( eval( lhs ) * eval( orhs_ ) );
            osres_  -= declsym( eval( lhs ) * eval( orhs_ ) );
            refres_ -= declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declsym( eval( olhs ) * eval( rhs_ ) );
            odres_  -= declsym( eval( olhs ) * eval( rhs_ ) );
            sres_   -= declsym( eval( olhs ) * eval( rhs_ ) );
            osres_  -= declsym( eval( olhs ) * eval( rhs_ ) );
            refres_ -= declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declsym( eval( olhs ) * eval( orhs_ ) );
            odres_  -= declsym( eval( olhs ) * eval( orhs_ ) );
            sres_   -= declsym( eval( olhs ) * eval( orhs_ ) );
            osres_  -= declsym( eval( olhs ) * eval( orhs_ ) );
            refres_ -= declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declsym multiplication with Schur product assignment
      //=====================================================================================

      // Declsym multiplication with Schur product assignment with the given matrices
      {
         test_  = "Declsym multiplication with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declsym( lhs * rhs_ );
            odres_  %= declsym( lhs * rhs_ );
            sres_   %= declsym( lhs * rhs_ );
            osres_  %= declsym( lhs * rhs_ );
            refres_ %= declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declsym( lhs * orhs_ );
            odres_  %= declsym( lhs * orhs_ );
            sres_   %= declsym( lhs * orhs_ );
            osres_  %= declsym( lhs * orhs_ );
            refres_ %= declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declsym( olhs * rhs_ );
            odres_  %= declsym( olhs * rhs_ );
            sres_   %= declsym( olhs * rhs_ );
            osres_  %= declsym( olhs * rhs_ );
            refres_ %= declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declsym( olhs * orhs_ );
            odres_  %= declsym( olhs * orhs_ );
            sres_   %= declsym( olhs * orhs_ );
            osres_  %= declsym( olhs * orhs_ );
            refres_ %= declsym( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declsym multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Declsym multiplication with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declsym( eval( lhs ) * eval( rhs_ ) );
            odres_  %= declsym( eval( lhs ) * eval( rhs_ ) );
            sres_   %= declsym( eval( lhs ) * eval( rhs_ ) );
            osres_  %= declsym( eval( lhs ) * eval( rhs_ ) );
            refres_ %= declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declsym( eval( lhs ) * eval( orhs_ ) );
            odres_  %= declsym( eval( lhs ) * eval( orhs_ ) );
            sres_   %= declsym( eval( lhs ) * eval( orhs_ ) );
            osres_  %= declsym( eval( lhs ) * eval( orhs_ ) );
            refres_ %= declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declsym( eval( olhs ) * eval( rhs_ ) );
            odres_  %= declsym( eval( olhs ) * eval( rhs_ ) );
            sres_   %= declsym( eval( olhs ) * eval( rhs_ ) );
            osres_  %= declsym( eval( olhs ) * eval( rhs_ ) );
            refres_ %= declsym( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declsym( eval( olhs ) * eval( orhs_ ) );
            odres_  %= declsym( eval( olhs ) * eval( orhs_ ) );
            sres_   %= declsym( eval( olhs ) * eval( orhs_ ) );
            osres_  %= declsym( eval( olhs ) * eval( orhs_ ) );
            refres_ %= declsym( eval( reflhs ) * eval( refrhs_ ) );
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
/*!\brief Testing the Hermitian sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the Hermitian matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment in the case that the
// left-hand side element type is the lower data type. In case any error resulting from
// the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testDeclHermOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION > 1 )
   {
      if( lhs_.rows() != rhs_.columns() )
         return;


      //=====================================================================================
      // Test-specific setup
      //=====================================================================================

      MT2  rhs;
      OMT2 orhs;
      RT2  refrhs;

      try {
         rhs    = ctrans( lhs_ );
         orhs   = ctrans( lhs_ );
         refrhs = ctrans( lhs_ );
      }
      catch( std::exception& ) {
         return;
      }


      //=====================================================================================
      // Declherm multiplication
      //=====================================================================================

      // Declherm multiplication with the given matrices
      {
         test_  = "Declherm multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = declherm( lhs_ * rhs );
            odres_  = declherm( lhs_ * rhs );
            sres_   = declherm( lhs_ * rhs );
            osres_  = declherm( lhs_ * rhs );
            refres_ = declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declherm( lhs_ * orhs );
            odres_  = declherm( lhs_ * orhs );
            sres_   = declherm( lhs_ * orhs );
            osres_  = declherm( lhs_ * orhs );
            refres_ = declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declherm( olhs_ * rhs );
            odres_  = declherm( olhs_ * rhs );
            sres_   = declherm( olhs_ * rhs );
            osres_  = declherm( olhs_ * rhs );
            refres_ = declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declherm( olhs_ * orhs );
            odres_  = declherm( olhs_ * orhs );
            sres_   = declherm( olhs_ * orhs );
            osres_  = declherm( olhs_ * orhs );
            refres_ = declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declherm multiplication with evaluated matrices
      {
         test_  = "Declherm multiplication with evaluated left-hand side matrix";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = declherm( eval( lhs_ ) * eval( rhs ) );
            odres_  = declherm( eval( lhs_ ) * eval( rhs ) );
            sres_   = declherm( eval( lhs_ ) * eval( rhs ) );
            osres_  = declherm( eval( lhs_ ) * eval( rhs ) );
            refres_ = declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declherm( eval( lhs_ ) * eval( orhs ) );
            odres_  = declherm( eval( lhs_ ) * eval( orhs ) );
            sres_   = declherm( eval( lhs_ ) * eval( orhs ) );
            osres_  = declherm( eval( lhs_ ) * eval( orhs ) );
            refres_ = declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declherm( eval( olhs_ ) * eval( rhs ) );
            odres_  = declherm( eval( olhs_ ) * eval( rhs ) );
            sres_   = declherm( eval( olhs_ ) * eval( rhs ) );
            osres_  = declherm( eval( olhs_ ) * eval( rhs ) );
            refres_ = declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declherm( eval( olhs_ ) * eval( orhs ) );
            odres_  = declherm( eval( olhs_ ) * eval( orhs ) );
            sres_   = declherm( eval( olhs_ ) * eval( orhs ) );
            osres_  = declherm( eval( olhs_ ) * eval( orhs ) );
            refres_ = declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declherm multiplication with addition assignment
      //=====================================================================================

      // Declherm multiplication with addition assignment with the given matrices
      {
         test_  = "Declherm multiplication with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declherm( lhs_ * rhs );
            odres_  += declherm( lhs_ * rhs );
            sres_   += declherm( lhs_ * rhs );
            osres_  += declherm( lhs_ * rhs );
            refres_ += declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declherm( lhs_ * orhs );
            odres_  += declherm( lhs_ * orhs );
            sres_   += declherm( lhs_ * orhs );
            osres_  += declherm( lhs_ * orhs );
            refres_ += declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declherm( olhs_ * rhs );
            odres_  += declherm( olhs_ * rhs );
            sres_   += declherm( olhs_ * rhs );
            osres_  += declherm( olhs_ * rhs );
            refres_ += declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declherm( olhs_ * orhs );
            odres_  += declherm( olhs_ * orhs );
            sres_   += declherm( olhs_ * orhs );
            osres_  += declherm( olhs_ * orhs );
            refres_ += declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declherm multiplication with addition assignment with evaluated matrices
      {
         test_  = "Declherm multiplication with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declherm( eval( lhs_ ) * eval( rhs ) );
            odres_  += declherm( eval( lhs_ ) * eval( rhs ) );
            sres_   += declherm( eval( lhs_ ) * eval( rhs ) );
            osres_  += declherm( eval( lhs_ ) * eval( rhs ) );
            refres_ += declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declherm( eval( lhs_ ) * eval( orhs ) );
            odres_  += declherm( eval( lhs_ ) * eval( orhs ) );
            sres_   += declherm( eval( lhs_ ) * eval( orhs ) );
            osres_  += declherm( eval( lhs_ ) * eval( orhs ) );
            refres_ += declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declherm( eval( olhs_ ) * eval( rhs ) );
            odres_  += declherm( eval( olhs_ ) * eval( rhs ) );
            sres_   += declherm( eval( olhs_ ) * eval( rhs ) );
            osres_  += declherm( eval( olhs_ ) * eval( rhs ) );
            refres_ += declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declherm( eval( olhs_ ) * eval( orhs ) );
            odres_  += declherm( eval( olhs_ ) * eval( orhs ) );
            sres_   += declherm( eval( olhs_ ) * eval( orhs ) );
            osres_  += declherm( eval( olhs_ ) * eval( orhs ) );
            refres_ += declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declherm multiplication with subtraction assignment
      //=====================================================================================

      // Declherm multiplication with subtraction assignment with the given matrices
      {
         test_  = "Declherm multiplication with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declherm( lhs_ * rhs );
            odres_  -= declherm( lhs_ * rhs );
            sres_   -= declherm( lhs_ * rhs );
            osres_  -= declherm( lhs_ * rhs );
            refres_ -= declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declherm( lhs_ * orhs );
            odres_  -= declherm( lhs_ * orhs );
            sres_   -= declherm( lhs_ * orhs );
            osres_  -= declherm( lhs_ * orhs );
            refres_ -= declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declherm( olhs_ * rhs );
            odres_  -= declherm( olhs_ * rhs );
            sres_   -= declherm( olhs_ * rhs );
            osres_  -= declherm( olhs_ * rhs );
            refres_ -= declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declherm( olhs_ * orhs );
            odres_  -= declherm( olhs_ * orhs );
            sres_   -= declherm( olhs_ * orhs );
            osres_  -= declherm( olhs_ * orhs );
            refres_ -= declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declherm multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Declherm multiplication with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declherm( eval( lhs_ ) * eval( rhs ) );
            odres_  -= declherm( eval( lhs_ ) * eval( rhs ) );
            sres_   -= declherm( eval( lhs_ ) * eval( rhs ) );
            osres_  -= declherm( eval( lhs_ ) * eval( rhs ) );
            refres_ -= declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declherm( eval( lhs_ ) * eval( orhs ) );
            odres_  -= declherm( eval( lhs_ ) * eval( orhs ) );
            sres_   -= declherm( eval( lhs_ ) * eval( orhs ) );
            osres_  -= declherm( eval( lhs_ ) * eval( orhs ) );
            refres_ -= declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declherm( eval( olhs_ ) * eval( rhs ) );
            odres_  -= declherm( eval( olhs_ ) * eval( rhs ) );
            sres_   -= declherm( eval( olhs_ ) * eval( rhs ) );
            osres_  -= declherm( eval( olhs_ ) * eval( rhs ) );
            refres_ -= declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declherm( eval( olhs_ ) * eval( orhs ) );
            odres_  -= declherm( eval( olhs_ ) * eval( orhs ) );
            sres_   -= declherm( eval( olhs_ ) * eval( orhs ) );
            osres_  -= declherm( eval( olhs_ ) * eval( orhs ) );
            refres_ -= declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declherm multiplication with Schur product assignment
      //=====================================================================================

      // Declherm multiplication with Schur product assignment with the given matrices
      {
         test_  = "Declherm multiplication with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declherm( lhs_ * rhs );
            odres_  %= declherm( lhs_ * rhs );
            sres_   %= declherm( lhs_ * rhs );
            osres_  %= declherm( lhs_ * rhs );
            refres_ %= declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declherm( lhs_ * orhs );
            odres_  %= declherm( lhs_ * orhs );
            sres_   %= declherm( lhs_ * orhs );
            osres_  %= declherm( lhs_ * orhs );
            refres_ %= declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declherm( olhs_ * rhs );
            odres_  %= declherm( olhs_ * rhs );
            sres_   %= declherm( olhs_ * rhs );
            osres_  %= declherm( olhs_ * rhs );
            refres_ %= declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declherm( olhs_ * orhs );
            odres_  %= declherm( olhs_ * orhs );
            sres_   %= declherm( olhs_ * orhs );
            osres_  %= declherm( olhs_ * orhs );
            refres_ %= declherm( reflhs_ * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declherm multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Declherm multiplication with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declherm( eval( lhs_ ) * eval( rhs ) );
            odres_  %= declherm( eval( lhs_ ) * eval( rhs ) );
            sres_   %= declherm( eval( lhs_ ) * eval( rhs ) );
            osres_  %= declherm( eval( lhs_ ) * eval( rhs ) );
            refres_ %= declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declherm( eval( lhs_ ) * eval( orhs ) );
            odres_  %= declherm( eval( lhs_ ) * eval( orhs ) );
            sres_   %= declherm( eval( lhs_ ) * eval( orhs ) );
            osres_  %= declherm( eval( lhs_ ) * eval( orhs ) );
            refres_ %= declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declherm( eval( olhs_ ) * eval( rhs ) );
            odres_  %= declherm( eval( olhs_ ) * eval( rhs ) );
            sres_   %= declherm( eval( olhs_ ) * eval( rhs ) );
            osres_  %= declherm( eval( olhs_ ) * eval( rhs ) );
            refres_ %= declherm( eval( reflhs_ ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declherm( eval( olhs_ ) * eval( orhs ) );
            odres_  %= declherm( eval( olhs_ ) * eval( orhs ) );
            sres_   %= declherm( eval( olhs_ ) * eval( orhs ) );
            osres_  %= declherm( eval( olhs_ ) * eval( orhs ) );
            refres_ %= declherm( eval( reflhs_ ) * eval( refrhs ) );
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
/*!\brief Testing the Hermitian sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the Hermitian matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment in the case that the
// right-hand side element type is the lower data type. In case any error resulting from
// the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testDeclHermOperation( blaze::FalseType )
{
#if BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION > 1 )
   {
      if( lhs_.rows() != rhs_.columns() )
         return;


      //=====================================================================================
      // Test-specific setup
      //=====================================================================================

      MT1  lhs;
      OMT1 olhs;
      RT1  reflhs;

      try {
         lhs    = ctrans( rhs_ );
         olhs   = ctrans( rhs_ );
         reflhs = ctrans( rhs_ );
      }
      catch( std::exception& ) {
         return;
      }


      //=====================================================================================
      // Declherm multiplication
      //=====================================================================================

      // Declherm multiplication with the given matrices
      {
         test_  = "Declherm multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = declherm( lhs * rhs_ );
            odres_  = declherm( lhs * rhs_ );
            sres_   = declherm( lhs * rhs_ );
            osres_  = declherm( lhs * rhs_ );
            refres_ = declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declherm( lhs * orhs_ );
            odres_  = declherm( lhs * orhs_ );
            sres_   = declherm( lhs * orhs_ );
            osres_  = declherm( lhs * orhs_ );
            refres_ = declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declherm( olhs * rhs_ );
            odres_  = declherm( olhs * rhs_ );
            sres_   = declherm( olhs * rhs_ );
            osres_  = declherm( olhs * rhs_ );
            refres_ = declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declherm( olhs * orhs_ );
            odres_  = declherm( olhs * orhs_ );
            sres_   = declherm( olhs * orhs_ );
            osres_  = declherm( olhs * orhs_ );
            refres_ = declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declherm multiplication with evaluated matrices
      {
         test_  = "Declherm multiplication with evaluated left-hand side matrix";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = declherm( eval( lhs ) * eval( rhs_ ) );
            odres_  = declherm( eval( lhs ) * eval( rhs_ ) );
            sres_   = declherm( eval( lhs ) * eval( rhs_ ) );
            osres_  = declherm( eval( lhs ) * eval( rhs_ ) );
            refres_ = declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declherm( eval( lhs ) * eval( orhs_ ) );
            odres_  = declherm( eval( lhs ) * eval( orhs_ ) );
            sres_   = declherm( eval( lhs ) * eval( orhs_ ) );
            osres_  = declherm( eval( lhs ) * eval( orhs_ ) );
            refres_ = declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declherm( eval( olhs ) * eval( rhs_ ) );
            odres_  = declherm( eval( olhs ) * eval( rhs_ ) );
            sres_   = declherm( eval( olhs ) * eval( rhs_ ) );
            osres_  = declherm( eval( olhs ) * eval( rhs_ ) );
            refres_ = declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declherm( eval( olhs ) * eval( orhs_ ) );
            odres_  = declherm( eval( olhs ) * eval( orhs_ ) );
            sres_   = declherm( eval( olhs ) * eval( orhs_ ) );
            osres_  = declherm( eval( olhs ) * eval( orhs_ ) );
            refres_ = declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declherm multiplication with addition assignment
      //=====================================================================================

      // Declherm multiplication with addition assignment with the given matrices
      {
         test_  = "Declherm multiplication with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declherm( lhs * rhs_ );
            odres_  += declherm( lhs * rhs_ );
            sres_   += declherm( lhs * rhs_ );
            osres_  += declherm( lhs * rhs_ );
            refres_ += declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declherm( lhs * orhs_ );
            odres_  += declherm( lhs * orhs_ );
            sres_   += declherm( lhs * orhs_ );
            osres_  += declherm( lhs * orhs_ );
            refres_ += declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declherm( olhs * rhs_ );
            odres_  += declherm( olhs * rhs_ );
            sres_   += declherm( olhs * rhs_ );
            osres_  += declherm( olhs * rhs_ );
            refres_ += declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declherm( olhs * orhs_ );
            odres_  += declherm( olhs * orhs_ );
            sres_   += declherm( olhs * orhs_ );
            osres_  += declherm( olhs * orhs_ );
            refres_ += declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declherm multiplication with addition assignment with evaluated matrices
      {
         test_  = "Declherm multiplication with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declherm( eval( lhs ) * eval( rhs_ ) );
            odres_  += declherm( eval( lhs ) * eval( rhs_ ) );
            sres_   += declherm( eval( lhs ) * eval( rhs_ ) );
            osres_  += declherm( eval( lhs ) * eval( rhs_ ) );
            refres_ += declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declherm( eval( lhs ) * eval( orhs_ ) );
            odres_  += declherm( eval( lhs ) * eval( orhs_ ) );
            sres_   += declherm( eval( lhs ) * eval( orhs_ ) );
            osres_  += declherm( eval( lhs ) * eval( orhs_ ) );
            refres_ += declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declherm( eval( olhs ) * eval( rhs_ ) );
            odres_  += declherm( eval( olhs ) * eval( rhs_ ) );
            sres_   += declherm( eval( olhs ) * eval( rhs_ ) );
            osres_  += declherm( eval( olhs ) * eval( rhs_ ) );
            refres_ += declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declherm( eval( olhs ) * eval( orhs_ ) );
            odres_  += declherm( eval( olhs ) * eval( orhs_ ) );
            sres_   += declherm( eval( olhs ) * eval( orhs_ ) );
            osres_  += declherm( eval( olhs ) * eval( orhs_ ) );
            refres_ += declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declherm multiplication with subtraction assignment
      //=====================================================================================

      // Declherm multiplication with subtraction assignment with the given matrices
      {
         test_  = "Declherm multiplication with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declherm( lhs * rhs_ );
            odres_  -= declherm( lhs * rhs_ );
            sres_   -= declherm( lhs * rhs_ );
            osres_  -= declherm( lhs * rhs_ );
            refres_ -= declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declherm( lhs * orhs_ );
            odres_  -= declherm( lhs * orhs_ );
            sres_   -= declherm( lhs * orhs_ );
            osres_  -= declherm( lhs * orhs_ );
            refres_ -= declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declherm( olhs * rhs_ );
            odres_  -= declherm( olhs * rhs_ );
            sres_   -= declherm( olhs * rhs_ );
            osres_  -= declherm( olhs * rhs_ );
            refres_ -= declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declherm( olhs * orhs_ );
            odres_  -= declherm( olhs * orhs_ );
            sres_   -= declherm( olhs * orhs_ );
            osres_  -= declherm( olhs * orhs_ );
            refres_ -= declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declherm multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Declherm multiplication with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declherm( eval( lhs ) * eval( rhs_ ) );
            odres_  -= declherm( eval( lhs ) * eval( rhs_ ) );
            sres_   -= declherm( eval( lhs ) * eval( rhs_ ) );
            osres_  -= declherm( eval( lhs ) * eval( rhs_ ) );
            refres_ -= declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declherm( eval( lhs ) * eval( orhs_ ) );
            odres_  -= declherm( eval( lhs ) * eval( orhs_ ) );
            sres_   -= declherm( eval( lhs ) * eval( orhs_ ) );
            osres_  -= declherm( eval( lhs ) * eval( orhs_ ) );
            refres_ -= declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declherm( eval( olhs ) * eval( rhs_ ) );
            odres_  -= declherm( eval( olhs ) * eval( rhs_ ) );
            sres_   -= declherm( eval( olhs ) * eval( rhs_ ) );
            osres_  -= declherm( eval( olhs ) * eval( rhs_ ) );
            refres_ -= declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declherm( eval( olhs ) * eval( orhs_ ) );
            odres_  -= declherm( eval( olhs ) * eval( orhs_ ) );
            sres_   -= declherm( eval( olhs ) * eval( orhs_ ) );
            osres_  -= declherm( eval( olhs ) * eval( orhs_ ) );
            refres_ -= declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declherm multiplication with Schur product assignment
      //=====================================================================================

      // Declherm multiplication with Schur product assignment with the given matrices
      {
         test_  = "Declherm multiplication with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declherm( lhs * rhs_ );
            odres_  %= declherm( lhs * rhs_ );
            sres_   %= declherm( lhs * rhs_ );
            osres_  %= declherm( lhs * rhs_ );
            refres_ %= declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declherm( lhs * orhs_ );
            odres_  %= declherm( lhs * orhs_ );
            sres_   %= declherm( lhs * orhs_ );
            osres_  %= declherm( lhs * orhs_ );
            refres_ %= declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declherm( olhs * rhs_ );
            odres_  %= declherm( olhs * rhs_ );
            sres_   %= declherm( olhs * rhs_ );
            osres_  %= declherm( olhs * rhs_ );
            refres_ %= declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declherm( olhs * orhs_ );
            odres_  %= declherm( olhs * orhs_ );
            sres_   %= declherm( olhs * orhs_ );
            osres_  %= declherm( olhs * orhs_ );
            refres_ %= declherm( reflhs * refrhs_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declherm multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Declherm multiplication with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declherm( eval( lhs ) * eval( rhs_ ) );
            odres_  %= declherm( eval( lhs ) * eval( rhs_ ) );
            sres_   %= declherm( eval( lhs ) * eval( rhs_ ) );
            osres_  %= declherm( eval( lhs ) * eval( rhs_ ) );
            refres_ %= declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declherm( eval( lhs ) * eval( orhs_ ) );
            odres_  %= declherm( eval( lhs ) * eval( orhs_ ) );
            sres_   %= declherm( eval( lhs ) * eval( orhs_ ) );
            osres_  %= declherm( eval( lhs ) * eval( orhs_ ) );
            refres_ %= declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declherm( eval( olhs ) * eval( rhs_ ) );
            odres_  %= declherm( eval( olhs ) * eval( rhs_ ) );
            sres_   %= declherm( eval( olhs ) * eval( rhs_ ) );
            osres_  %= declherm( eval( olhs ) * eval( rhs_ ) );
            refres_ %= declherm( eval( reflhs ) * eval( refrhs_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declherm( eval( olhs ) * eval( orhs_ ) );
            odres_  %= declherm( eval( olhs ) * eval( orhs_ ) );
            sres_   %= declherm( eval( olhs ) * eval( orhs_ ) );
            osres_  %= declherm( eval( olhs ) * eval( orhs_ ) );
            refres_ %= declherm( eval( reflhs ) * eval( refrhs_ ) );
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
/*!\brief Testing the lower sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the lower matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testDeclLowOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLLOW_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLLOW_OPERATION > 1 )
   {
      if( lhs_.rows() != rhs_.columns() )
         return;


      //=====================================================================================
      // Test-specific setup of the left-hand side operand
      //=====================================================================================

      MT1 lhs( lhs_ );

      blaze::resetUpper( lhs );

      OMT1 olhs  ( lhs );
      RT1  reflhs( lhs );


      //=====================================================================================
      // Test-specific setup of the right-hand side operand
      //=====================================================================================

      MT2 rhs( rhs_ );

      blaze::resetUpper( rhs );

      OMT2 orhs  ( rhs );
      RT2  refrhs( rhs );


      //=====================================================================================
      // Decllow multiplication
      //=====================================================================================

      // Decllow multiplication with the given matrices
      {
         test_  = "Decllow multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = decllow( lhs * rhs );
            odres_  = decllow( lhs * rhs );
            sres_   = decllow( lhs * rhs );
            osres_  = decllow( lhs * rhs );
            refres_ = decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = decllow( lhs * orhs );
            odres_  = decllow( lhs * orhs );
            sres_   = decllow( lhs * orhs );
            osres_  = decllow( lhs * orhs );
            refres_ = decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = decllow( olhs * rhs );
            odres_  = decllow( olhs * rhs );
            sres_   = decllow( olhs * rhs );
            osres_  = decllow( olhs * rhs );
            refres_ = decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = decllow( olhs * orhs );
            odres_  = decllow( olhs * orhs );
            sres_   = decllow( olhs * orhs );
            osres_  = decllow( olhs * orhs );
            refres_ = decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decllow multiplication with evaluated matrices
      {
         test_  = "Decllow multiplication with evaluated left-hand side matrix";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = decllow( eval( lhs ) * eval( rhs ) );
            odres_  = decllow( eval( lhs ) * eval( rhs ) );
            sres_   = decllow( eval( lhs ) * eval( rhs ) );
            osres_  = decllow( eval( lhs ) * eval( rhs ) );
            refres_ = decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = decllow( eval( lhs ) * eval( orhs ) );
            odres_  = decllow( eval( lhs ) * eval( orhs ) );
            sres_   = decllow( eval( lhs ) * eval( orhs ) );
            osres_  = decllow( eval( lhs ) * eval( orhs ) );
            refres_ = decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = decllow( eval( olhs ) * eval( rhs ) );
            odres_  = decllow( eval( olhs ) * eval( rhs ) );
            sres_   = decllow( eval( olhs ) * eval( rhs ) );
            osres_  = decllow( eval( olhs ) * eval( rhs ) );
            refres_ = decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = decllow( eval( olhs ) * eval( orhs ) );
            odres_  = decllow( eval( olhs ) * eval( orhs ) );
            sres_   = decllow( eval( olhs ) * eval( orhs ) );
            osres_  = decllow( eval( olhs ) * eval( orhs ) );
            refres_ = decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Decllow multiplication with addition assignment
      //=====================================================================================

      // Decllow multiplication with addition assignment with the given matrices
      {
         test_  = "Decllow multiplication with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += decllow( lhs * rhs );
            odres_  += decllow( lhs * rhs );
            sres_   += decllow( lhs * rhs );
            osres_  += decllow( lhs * rhs );
            refres_ += decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += decllow( lhs * orhs );
            odres_  += decllow( lhs * orhs );
            sres_   += decllow( lhs * orhs );
            osres_  += decllow( lhs * orhs );
            refres_ += decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += decllow( olhs * rhs );
            odres_  += decllow( olhs * rhs );
            sres_   += decllow( olhs * rhs );
            osres_  += decllow( olhs * rhs );
            refres_ += decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += decllow( olhs * orhs );
            odres_  += decllow( olhs * orhs );
            sres_   += decllow( olhs * orhs );
            osres_  += decllow( olhs * orhs );
            refres_ += decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decllow multiplication with addition assignment with evaluated matrices
      {
         test_  = "Decllow multiplication with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += decllow( eval( lhs ) * eval( rhs ) );
            odres_  += decllow( eval( lhs ) * eval( rhs ) );
            sres_   += decllow( eval( lhs ) * eval( rhs ) );
            osres_  += decllow( eval( lhs ) * eval( rhs ) );
            refres_ += decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += decllow( eval( lhs ) * eval( orhs ) );
            odres_  += decllow( eval( lhs ) * eval( orhs ) );
            sres_   += decllow( eval( lhs ) * eval( orhs ) );
            osres_  += decllow( eval( lhs ) * eval( orhs ) );
            refres_ += decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += decllow( eval( olhs ) * eval( rhs ) );
            odres_  += decllow( eval( olhs ) * eval( rhs ) );
            sres_   += decllow( eval( olhs ) * eval( rhs ) );
            osres_  += decllow( eval( olhs ) * eval( rhs ) );
            refres_ += decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += decllow( eval( olhs ) * eval( orhs ) );
            odres_  += decllow( eval( olhs ) * eval( orhs ) );
            sres_   += decllow( eval( olhs ) * eval( orhs ) );
            osres_  += decllow( eval( olhs ) * eval( orhs ) );
            refres_ += decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Decllow multiplication with subtraction assignment
      //=====================================================================================

      // Decllow multiplication with subtraction assignment with the given matrices
      {
         test_  = "Decllow multiplication with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= decllow( lhs * rhs );
            odres_  -= decllow( lhs * rhs );
            sres_   -= decllow( lhs * rhs );
            osres_  -= decllow( lhs * rhs );
            refres_ -= decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= decllow( lhs * orhs );
            odres_  -= decllow( lhs * orhs );
            sres_   -= decllow( lhs * orhs );
            osres_  -= decllow( lhs * orhs );
            refres_ -= decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= decllow( olhs * rhs );
            odres_  -= decllow( olhs * rhs );
            sres_   -= decllow( olhs * rhs );
            osres_  -= decllow( olhs * rhs );
            refres_ -= decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= decllow( olhs * orhs );
            odres_  -= decllow( olhs * orhs );
            sres_   -= decllow( olhs * orhs );
            osres_  -= decllow( olhs * orhs );
            refres_ -= decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decllow multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Decllow multiplication with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= decllow( eval( lhs ) * eval( rhs ) );
            odres_  -= decllow( eval( lhs ) * eval( rhs ) );
            sres_   -= decllow( eval( lhs ) * eval( rhs ) );
            osres_  -= decllow( eval( lhs ) * eval( rhs ) );
            refres_ -= decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= decllow( eval( lhs ) * eval( orhs ) );
            odres_  -= decllow( eval( lhs ) * eval( orhs ) );
            sres_   -= decllow( eval( lhs ) * eval( orhs ) );
            osres_  -= decllow( eval( lhs ) * eval( orhs ) );
            refres_ -= decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= decllow( eval( olhs ) * eval( rhs ) );
            odres_  -= decllow( eval( olhs ) * eval( rhs ) );
            sres_   -= decllow( eval( olhs ) * eval( rhs ) );
            osres_  -= decllow( eval( olhs ) * eval( rhs ) );
            refres_ -= decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= decllow( eval( olhs ) * eval( orhs ) );
            odres_  -= decllow( eval( olhs ) * eval( orhs ) );
            sres_   -= decllow( eval( olhs ) * eval( orhs ) );
            osres_  -= decllow( eval( olhs ) * eval( orhs ) );
            refres_ -= decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Decllow multiplication with Schur product assignment
      //=====================================================================================

      // Decllow multiplication with Schur product assignment with the given matrices
      {
         test_  = "Decllow multiplication with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= decllow( lhs * rhs );
            odres_  %= decllow( lhs * rhs );
            sres_   %= decllow( lhs * rhs );
            osres_  %= decllow( lhs * rhs );
            refres_ %= decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= decllow( lhs * orhs );
            odres_  %= decllow( lhs * orhs );
            sres_   %= decllow( lhs * orhs );
            osres_  %= decllow( lhs * orhs );
            refres_ %= decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= decllow( olhs * rhs );
            odres_  %= decllow( olhs * rhs );
            sres_   %= decllow( olhs * rhs );
            osres_  %= decllow( olhs * rhs );
            refres_ %= decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= decllow( olhs * orhs );
            odres_  %= decllow( olhs * orhs );
            sres_   %= decllow( olhs * orhs );
            osres_  %= decllow( olhs * orhs );
            refres_ %= decllow( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decllow multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Decllow multiplication with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= decllow( eval( lhs ) * eval( rhs ) );
            odres_  %= decllow( eval( lhs ) * eval( rhs ) );
            sres_   %= decllow( eval( lhs ) * eval( rhs ) );
            osres_  %= decllow( eval( lhs ) * eval( rhs ) );
            refres_ %= decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= decllow( eval( lhs ) * eval( orhs ) );
            odres_  %= decllow( eval( lhs ) * eval( orhs ) );
            sres_   %= decllow( eval( lhs ) * eval( orhs ) );
            osres_  %= decllow( eval( lhs ) * eval( orhs ) );
            refres_ %= decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= decllow( eval( olhs ) * eval( rhs ) );
            odres_  %= decllow( eval( olhs ) * eval( rhs ) );
            sres_   %= decllow( eval( olhs ) * eval( rhs ) );
            osres_  %= decllow( eval( olhs ) * eval( rhs ) );
            refres_ %= decllow( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= decllow( eval( olhs ) * eval( orhs ) );
            odres_  %= decllow( eval( olhs ) * eval( orhs ) );
            sres_   %= decllow( eval( olhs ) * eval( orhs ) );
            osres_  %= decllow( eval( olhs ) * eval( orhs ) );
            refres_ %= decllow( eval( reflhs ) * eval( refrhs ) );
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
/*!\brief Skipping the lower sparse matrix/sparse matrix multiplication.
//
// \return void
//
// This function is called in case the lower matrix/matrix multiplication operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testDeclLowOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the upper sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the upper matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testDeclUppOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLUPP_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLUPP_OPERATION > 1 )
   {
      if( lhs_.rows() != rhs_.columns() )
         return;


      //=====================================================================================
      // Test-specific setup of the left-hand side operand
      //=====================================================================================

      MT1 lhs( lhs_ );

      blaze::resetLower( lhs );

      OMT1 olhs  ( lhs );
      RT1  reflhs( lhs );


      //=====================================================================================
      // Test-specific setup of the right-hand side operand
      //=====================================================================================

      MT2 rhs( rhs_ );

      blaze::resetLower( rhs );

      OMT2 orhs  ( rhs );
      RT2  refrhs( rhs );


      //=====================================================================================
      // Declupp multiplication
      //=====================================================================================

      // Declupp multiplication with the given matrices
      {
         test_  = "Declupp multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = declupp( lhs * rhs );
            odres_  = declupp( lhs * rhs );
            sres_   = declupp( lhs * rhs );
            osres_  = declupp( lhs * rhs );
            refres_ = declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declupp( lhs * orhs );
            odres_  = declupp( lhs * orhs );
            sres_   = declupp( lhs * orhs );
            osres_  = declupp( lhs * orhs );
            refres_ = declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declupp( olhs * rhs );
            odres_  = declupp( olhs * rhs );
            sres_   = declupp( olhs * rhs );
            osres_  = declupp( olhs * rhs );
            refres_ = declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declupp( olhs * orhs );
            odres_  = declupp( olhs * orhs );
            sres_   = declupp( olhs * orhs );
            osres_  = declupp( olhs * orhs );
            refres_ = declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declupp multiplication with evaluated matrices
      {
         test_  = "Declupp multiplication with evaluated left-hand side matrix";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = declupp( eval( lhs ) * eval( rhs ) );
            odres_  = declupp( eval( lhs ) * eval( rhs ) );
            sres_   = declupp( eval( lhs ) * eval( rhs ) );
            osres_  = declupp( eval( lhs ) * eval( rhs ) );
            refres_ = declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = declupp( eval( lhs ) * eval( orhs ) );
            odres_  = declupp( eval( lhs ) * eval( orhs ) );
            sres_   = declupp( eval( lhs ) * eval( orhs ) );
            osres_  = declupp( eval( lhs ) * eval( orhs ) );
            refres_ = declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = declupp( eval( olhs ) * eval( rhs ) );
            odres_  = declupp( eval( olhs ) * eval( rhs ) );
            sres_   = declupp( eval( olhs ) * eval( rhs ) );
            osres_  = declupp( eval( olhs ) * eval( rhs ) );
            refres_ = declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = declupp( eval( olhs ) * eval( orhs ) );
            odres_  = declupp( eval( olhs ) * eval( orhs ) );
            sres_   = declupp( eval( olhs ) * eval( orhs ) );
            osres_  = declupp( eval( olhs ) * eval( orhs ) );
            refres_ = declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declupp multiplication with addition assignment
      //=====================================================================================

      // Declupp multiplication with addition assignment with the given matrices
      {
         test_  = "Declupp multiplication with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declupp( lhs * rhs );
            odres_  += declupp( lhs * rhs );
            sres_   += declupp( lhs * rhs );
            osres_  += declupp( lhs * rhs );
            refres_ += declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declupp( lhs * orhs );
            odres_  += declupp( lhs * orhs );
            sres_   += declupp( lhs * orhs );
            osres_  += declupp( lhs * orhs );
            refres_ += declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declupp( olhs * rhs );
            odres_  += declupp( olhs * rhs );
            sres_   += declupp( olhs * rhs );
            osres_  += declupp( olhs * rhs );
            refres_ += declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declupp( olhs * orhs );
            odres_  += declupp( olhs * orhs );
            sres_   += declupp( olhs * orhs );
            osres_  += declupp( olhs * orhs );
            refres_ += declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declupp multiplication with addition assignment with evaluated matrices
      {
         test_  = "Declupp multiplication with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += declupp( eval( lhs ) * eval( rhs ) );
            odres_  += declupp( eval( lhs ) * eval( rhs ) );
            sres_   += declupp( eval( lhs ) * eval( rhs ) );
            osres_  += declupp( eval( lhs ) * eval( rhs ) );
            refres_ += declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += declupp( eval( lhs ) * eval( orhs ) );
            odres_  += declupp( eval( lhs ) * eval( orhs ) );
            sres_   += declupp( eval( lhs ) * eval( orhs ) );
            osres_  += declupp( eval( lhs ) * eval( orhs ) );
            refres_ += declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += declupp( eval( olhs ) * eval( rhs ) );
            odres_  += declupp( eval( olhs ) * eval( rhs ) );
            sres_   += declupp( eval( olhs ) * eval( rhs ) );
            osres_  += declupp( eval( olhs ) * eval( rhs ) );
            refres_ += declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += declupp( eval( olhs ) * eval( orhs ) );
            odres_  += declupp( eval( olhs ) * eval( orhs ) );
            sres_   += declupp( eval( olhs ) * eval( orhs ) );
            osres_  += declupp( eval( olhs ) * eval( orhs ) );
            refres_ += declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declupp multiplication with subtraction assignment
      //=====================================================================================

      // Declupp multiplication with subtraction assignment with the given matrices
      {
         test_  = "Declupp multiplication with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declupp( lhs * rhs );
            odres_  -= declupp( lhs * rhs );
            sres_   -= declupp( lhs * rhs );
            osres_  -= declupp( lhs * rhs );
            refres_ -= declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declupp( lhs * orhs );
            odres_  -= declupp( lhs * orhs );
            sres_   -= declupp( lhs * orhs );
            osres_  -= declupp( lhs * orhs );
            refres_ -= declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declupp( olhs * rhs );
            odres_  -= declupp( olhs * rhs );
            sres_   -= declupp( olhs * rhs );
            osres_  -= declupp( olhs * rhs );
            refres_ -= declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declupp( olhs * orhs );
            odres_  -= declupp( olhs * orhs );
            sres_   -= declupp( olhs * orhs );
            osres_  -= declupp( olhs * orhs );
            refres_ -= declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declupp multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Declupp multiplication with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= declupp( eval( lhs ) * eval( rhs ) );
            odres_  -= declupp( eval( lhs ) * eval( rhs ) );
            sres_   -= declupp( eval( lhs ) * eval( rhs ) );
            osres_  -= declupp( eval( lhs ) * eval( rhs ) );
            refres_ -= declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= declupp( eval( lhs ) * eval( orhs ) );
            odres_  -= declupp( eval( lhs ) * eval( orhs ) );
            sres_   -= declupp( eval( lhs ) * eval( orhs ) );
            osres_  -= declupp( eval( lhs ) * eval( orhs ) );
            refres_ -= declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= declupp( eval( olhs ) * eval( rhs ) );
            odres_  -= declupp( eval( olhs ) * eval( rhs ) );
            sres_   -= declupp( eval( olhs ) * eval( rhs ) );
            osres_  -= declupp( eval( olhs ) * eval( rhs ) );
            refres_ -= declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= declupp( eval( olhs ) * eval( orhs ) );
            odres_  -= declupp( eval( olhs ) * eval( orhs ) );
            sres_   -= declupp( eval( olhs ) * eval( orhs ) );
            osres_  -= declupp( eval( olhs ) * eval( orhs ) );
            refres_ -= declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Declupp multiplication with Schur product assignment
      //=====================================================================================

      // Declupp multiplication with Schur product assignment with the given matrices
      {
         test_  = "Declupp multiplication with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declupp( lhs * rhs );
            odres_  %= declupp( lhs * rhs );
            sres_   %= declupp( lhs * rhs );
            osres_  %= declupp( lhs * rhs );
            refres_ %= declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declupp( lhs * orhs );
            odres_  %= declupp( lhs * orhs );
            sres_   %= declupp( lhs * orhs );
            osres_  %= declupp( lhs * orhs );
            refres_ %= declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declupp( olhs * rhs );
            odres_  %= declupp( olhs * rhs );
            sres_   %= declupp( olhs * rhs );
            osres_  %= declupp( olhs * rhs );
            refres_ %= declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declupp( olhs * orhs );
            odres_  %= declupp( olhs * orhs );
            sres_   %= declupp( olhs * orhs );
            osres_  %= declupp( olhs * orhs );
            refres_ %= declupp( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Declupp multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Declupp multiplication with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            dres_   %= declupp( eval( lhs ) * eval( rhs ) );
            odres_  %= declupp( eval( lhs ) * eval( rhs ) );
            sres_   %= declupp( eval( lhs ) * eval( rhs ) );
            osres_  %= declupp( eval( lhs ) * eval( rhs ) );
            refres_ %= declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= declupp( eval( lhs ) * eval( orhs ) );
            odres_  %= declupp( eval( lhs ) * eval( orhs ) );
            sres_   %= declupp( eval( lhs ) * eval( orhs ) );
            osres_  %= declupp( eval( lhs ) * eval( orhs ) );
            refres_ %= declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= declupp( eval( olhs ) * eval( rhs ) );
            odres_  %= declupp( eval( olhs ) * eval( rhs ) );
            sres_   %= declupp( eval( olhs ) * eval( rhs ) );
            osres_  %= declupp( eval( olhs ) * eval( rhs ) );
            refres_ %= declupp( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= declupp( eval( olhs ) * eval( orhs ) );
            odres_  %= declupp( eval( olhs ) * eval( orhs ) );
            sres_   %= declupp( eval( olhs ) * eval( orhs ) );
            osres_  %= declupp( eval( olhs ) * eval( orhs ) );
            refres_ %= declupp( eval( reflhs ) * eval( refrhs ) );
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
/*!\brief Skipping the upper sparse matrix/sparse matrix multiplication.
//
// \return void
//
// This function is called in case the upper matrix/matrix multiplication operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testDeclUppOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the diagonal sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the diagonal matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testDeclDiagOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_DECLDIAG_OPERATION
   if( BLAZETEST_MATHTEST_TEST_DECLDIAG_OPERATION > 1 )
   {
      if( lhs_.rows() != rhs_.columns() )
         return;


      //=====================================================================================
      // Test-specific setup of the left-hand side operand
      //=====================================================================================

      MT1 lhs( lhs_ );

      blaze::resetLower( lhs );
      blaze::resetUpper( lhs );

      OMT1 olhs  ( lhs );
      RT1  reflhs( lhs );


      //=====================================================================================
      // Test-specific setup of the right-hand side operand
      //=====================================================================================

      MT2 rhs( rhs_ );

      blaze::resetLower( rhs );
      blaze::resetUpper( rhs );

      OMT2 orhs  ( rhs );
      RT2  refrhs( rhs );


      //=====================================================================================
      // Decldiag multiplication
      //=====================================================================================

      // Decldiag multiplication with the given matrices
      {
         test_  = "Decldiag multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = decldiag( lhs * rhs );
            odres_  = decldiag( lhs * rhs );
            sres_   = decldiag( lhs * rhs );
            osres_  = decldiag( lhs * rhs );
            refres_ = decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = decldiag( lhs * orhs );
            odres_  = decldiag( lhs * orhs );
            sres_   = decldiag( lhs * orhs );
            osres_  = decldiag( lhs * orhs );
            refres_ = decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = decldiag( olhs * rhs );
            odres_  = decldiag( olhs * rhs );
            sres_   = decldiag( olhs * rhs );
            osres_  = decldiag( olhs * rhs );
            refres_ = decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = decldiag( olhs * orhs );
            odres_  = decldiag( olhs * orhs );
            sres_   = decldiag( olhs * orhs );
            osres_  = decldiag( olhs * orhs );
            refres_ = decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decldiag multiplication with evaluated matrices
      {
         test_  = "Decldiag multiplication with evaluated left-hand side matrix";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            dres_   = decldiag( eval( lhs ) * eval( rhs ) );
            odres_  = decldiag( eval( lhs ) * eval( rhs ) );
            sres_   = decldiag( eval( lhs ) * eval( rhs ) );
            osres_  = decldiag( eval( lhs ) * eval( rhs ) );
            refres_ = decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   = decldiag( eval( lhs ) * eval( orhs ) );
            odres_  = decldiag( eval( lhs ) * eval( orhs ) );
            sres_   = decldiag( eval( lhs ) * eval( orhs ) );
            osres_  = decldiag( eval( lhs ) * eval( orhs ) );
            refres_ = decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   = decldiag( eval( olhs ) * eval( rhs ) );
            odres_  = decldiag( eval( olhs ) * eval( rhs ) );
            sres_   = decldiag( eval( olhs ) * eval( rhs ) );
            osres_  = decldiag( eval( olhs ) * eval( rhs ) );
            refres_ = decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   = decldiag( eval( olhs ) * eval( orhs ) );
            odres_  = decldiag( eval( olhs ) * eval( orhs ) );
            sres_   = decldiag( eval( olhs ) * eval( orhs ) );
            osres_  = decldiag( eval( olhs ) * eval( orhs ) );
            refres_ = decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Decldiag multiplication with addition assignment
      //=====================================================================================

      // Decldiag multiplication with addition assignment with the given matrices
      {
         test_  = "Decldiag multiplication with addition assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += decldiag( lhs * orhs );
            odres_  += decldiag( lhs * orhs );
            sres_   += decldiag( lhs * orhs );
            osres_  += decldiag( lhs * orhs );
            refres_ += decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += decldiag( olhs * rhs );
            odres_  += decldiag( olhs * rhs );
            sres_   += decldiag( olhs * rhs );
            osres_  += decldiag( olhs * rhs );
            refres_ += decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += decldiag( olhs * orhs );
            odres_  += decldiag( olhs * orhs );
            sres_   += decldiag( olhs * orhs );
            osres_  += decldiag( olhs * orhs );
            refres_ += decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decldiag multiplication with addition assignment with evaluated matrices
      {
         test_  = "Decldiag multiplication with addition assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   += decldiag( eval( lhs ) * eval( orhs ) );
            odres_  += decldiag( eval( lhs ) * eval( orhs ) );
            sres_   += decldiag( eval( lhs ) * eval( orhs ) );
            osres_  += decldiag( eval( lhs ) * eval( orhs ) );
            refres_ += decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   += decldiag( eval( olhs ) * eval( rhs ) );
            odres_  += decldiag( eval( olhs ) * eval( rhs ) );
            sres_   += decldiag( eval( olhs ) * eval( rhs ) );
            osres_  += decldiag( eval( olhs ) * eval( rhs ) );
            refres_ += decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   += decldiag( eval( olhs ) * eval( orhs ) );
            odres_  += decldiag( eval( olhs ) * eval( orhs ) );
            sres_   += decldiag( eval( olhs ) * eval( orhs ) );
            osres_  += decldiag( eval( olhs ) * eval( orhs ) );
            refres_ += decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Decldiag multiplication with subtraction assignment
      //=====================================================================================

      // Decldiag multiplication with subtraction assignment with the given matrices
      {
         test_  = "Decldiag multiplication with subtraction assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= decldiag( lhs * orhs );
            odres_  -= decldiag( lhs * orhs );
            sres_   -= decldiag( lhs * orhs );
            osres_  -= decldiag( lhs * orhs );
            refres_ -= decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= decldiag( olhs * rhs );
            odres_  -= decldiag( olhs * rhs );
            sres_   -= decldiag( olhs * rhs );
            osres_  -= decldiag( olhs * rhs );
            refres_ -= decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= decldiag( olhs * orhs );
            odres_  -= decldiag( olhs * orhs );
            sres_   -= decldiag( olhs * orhs );
            osres_  -= decldiag( olhs * orhs );
            refres_ -= decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decldiag multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Decldiag multiplication with subtraction assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   -= decldiag( eval( lhs ) * eval( orhs ) );
            odres_  -= decldiag( eval( lhs ) * eval( orhs ) );
            sres_   -= decldiag( eval( lhs ) * eval( orhs ) );
            osres_  -= decldiag( eval( lhs ) * eval( orhs ) );
            refres_ -= decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   -= decldiag( eval( olhs ) * eval( rhs ) );
            odres_  -= decldiag( eval( olhs ) * eval( rhs ) );
            sres_   -= decldiag( eval( olhs ) * eval( rhs ) );
            osres_  -= decldiag( eval( olhs ) * eval( rhs ) );
            refres_ -= decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   -= decldiag( eval( olhs ) * eval( orhs ) );
            odres_  -= decldiag( eval( olhs ) * eval( orhs ) );
            sres_   -= decldiag( eval( olhs ) * eval( orhs ) );
            osres_  -= decldiag( eval( olhs ) * eval( orhs ) );
            refres_ -= decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Decldiag multiplication with Schur product assignment
      //=====================================================================================

      // Decldiag multiplication with Schur product assignment with the given matrices
      {
         test_  = "Decldiag multiplication with Schur product assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= decldiag( lhs * orhs );
            odres_  %= decldiag( lhs * orhs );
            sres_   %= decldiag( lhs * orhs );
            osres_  %= decldiag( lhs * orhs );
            refres_ %= decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= decldiag( olhs * rhs );
            odres_  %= decldiag( olhs * rhs );
            sres_   %= decldiag( olhs * rhs );
            osres_  %= decldiag( olhs * rhs );
            refres_ %= decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= decldiag( olhs * orhs );
            odres_  %= decldiag( olhs * orhs );
            sres_   %= decldiag( olhs * orhs );
            osres_  %= decldiag( olhs * orhs );
            refres_ %= decldiag( reflhs * refrhs );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Decldiag multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Decldiag multiplication with Schur product assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            dres_   %= decldiag( eval( lhs ) * eval( orhs ) );
            odres_  %= decldiag( eval( lhs ) * eval( orhs ) );
            sres_   %= decldiag( eval( lhs ) * eval( orhs ) );
            osres_  %= decldiag( eval( lhs ) * eval( orhs ) );
            refres_ %= decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            dres_   %= decldiag( eval( olhs ) * eval( rhs ) );
            odres_  %= decldiag( eval( olhs ) * eval( rhs ) );
            sres_   %= decldiag( eval( olhs ) * eval( rhs ) );
            osres_  %= decldiag( eval( olhs ) * eval( rhs ) );
            refres_ %= decldiag( eval( reflhs ) * eval( refrhs ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            dres_   %= decldiag( eval( olhs ) * eval( orhs ) );
            odres_  %= decldiag( eval( olhs ) * eval( orhs ) );
            sres_   %= decldiag( eval( olhs ) * eval( orhs ) );
            osres_  %= decldiag( eval( olhs ) * eval( orhs ) );
            refres_ %= decldiag( eval( reflhs ) * eval( refrhs ) );
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
/*!\brief Skipping the diagonal sparse matrix/sparse matrix multiplication.
//
// \return void
//
// This function is called in case the diagonal matrix/matrix multiplication operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testDeclDiagOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the submatrix-wise sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the submatrix-wise matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testSubmatrixOperation()
{
#if BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION > 1 )
   {
      if( lhs_.rows() == 0UL || rhs_.columns() == 0UL )
         return;


      //=====================================================================================
      // Submatrix-wise multiplication
      //=====================================================================================

      // Submatrix-wise multiplication with the given matrices
      {
         test_  = "Submatrix-wise multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( reflhs_ * refrhs_, row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) = submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( reflhs_ * refrhs_, row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) = submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( reflhs_ * refrhs_, row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) = submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( reflhs_ * refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Submatrix-wise multiplication with evaluated matrices
      {
         test_  = "Submatrix-wise multiplication with evaluated matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) = submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) = submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Submatrix-wise multiplication with addition assignment
      //=====================================================================================

      // Submatrix-wise multiplication with addition assignment with the given matrices
      {
         test_  = "Submatrix-wise multiplication with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( reflhs_ * refrhs_, row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) += submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( reflhs_ * refrhs_, row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) += submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( reflhs_ * refrhs_, row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) += submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( reflhs_ * refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Submatrix-wise multiplication with addition assignment with evaluated matrices
      {
         test_  = "Submatrix-wise multiplication with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) += submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) += submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Submatrix-wise multiplication with subtraction assignment
      //=====================================================================================

      // Submatrix-wise multiplication with subtraction assignment with the given matrices
      {
         test_  = "Submatrix-wise multiplication with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( reflhs_ * refrhs_, row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) -= submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( reflhs_ * refrhs_, row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) -= submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( reflhs_ * refrhs_, row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) -= submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( reflhs_ * refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Submatrix-wise multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Submatrix-wise multiplication with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) -= submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) -= submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Submatrix-wise multiplication with Schur product assignment
      //=====================================================================================

      // Submatrix-wise multiplication with Schur product assignment with the given matrices
      {
         test_  = "Submatrix-wise multiplication with Schur product assignment with the given matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( lhs_ * rhs_      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( reflhs_ * refrhs_, row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) %= submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( lhs_ * orhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( reflhs_ * refrhs_, row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) %= submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( olhs_ * rhs_     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( reflhs_ * refrhs_, row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) %= submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( olhs_ * orhs_    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( reflhs_ * refrhs_, row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Submatrix-wise multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Submatrix-wise multiplication with Schur product assignment with evaluated matrices";
         error_ = "Failed Schur product assignment operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<lhs_.rows(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, lhs_.rows() - row );
               for( size_t column=0UL, n=0UL; column<rhs_.columns(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, rhs_.columns() - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( eval( lhs_ ) * eval( rhs_ )      , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) %= submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( eval( lhs_ ) * eval( orhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) %= submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( eval( olhs_ ) * eval( rhs_ )     , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
                  submatrix( dres_  , row, column, m, n ) %= submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( eval( olhs_ ) * eval( orhs_ )    , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( eval( reflhs_ ) * eval( refrhs_ ), row, column, m, n );
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
            auto sm = submatrix( lhs_ * rhs_, 1UL, 0UL, lhs_.rows(), rhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( lhs_ * orhs_, 1UL, 0UL, lhs_.rows(), orhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( olhs_ * rhs_, 1UL, 0UL, olhs_.rows(), rhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( olhs_ * orhs_, 1UL, 0UL, olhs_.rows(), orhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
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
            auto sm = submatrix( lhs_ * rhs_, 0UL, 1UL, lhs_.rows(), rhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( lhs_ * orhs_, 0UL, 1UL, lhs_.rows(), orhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( olhs_ * rhs_, 0UL, 1UL, olhs_.rows(), rhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( olhs_ * orhs_, 0UL, 1UL, olhs_.rows(), orhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
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
            auto sm = submatrix( lhs_ * rhs_, lhs_.rows(), 0UL, 1UL, rhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( lhs_ * orhs_, lhs_.rows(), 0UL, 1UL, orhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( olhs_ * rhs_, olhs_.rows(), 0UL, 1UL, rhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( olhs_ * orhs_, olhs_.rows(), 0UL, 1UL, orhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
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
            auto sm = submatrix( lhs_ * rhs_, 0UL, rhs_.columns(), lhs_.rows(), 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( lhs_ * orhs_, 0UL, orhs_.columns(), lhs_.rows(), 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( olhs_ * rhs_, 0UL, rhs_.columns(), olhs_.rows(), 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( olhs_ * orhs_, 0UL, orhs_.columns(), olhs_.rows(), 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
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
/*!\brief Testing the row-wise sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the row-wise matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testRowOperation()
{
#if BLAZETEST_MATHTEST_TEST_ROW_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ROW_OPERATION > 1 )
   {
      if( lhs_.rows() == 0UL )
         return;


      //=====================================================================================
      // Row-wise multiplication
      //=====================================================================================

      // Row-wise multiplication with the given matrices
      {
         test_  = "Row-wise multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( lhs_ * rhs_, i );
               row( odres_ , i ) = row( lhs_ * rhs_, i );
               row( sres_  , i ) = row( lhs_ * rhs_, i );
               row( osres_ , i ) = row( lhs_ * rhs_, i );
               row( refres_, i ) = row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( lhs_ * orhs_, i );
               row( odres_ , i ) = row( lhs_ * orhs_, i );
               row( sres_  , i ) = row( lhs_ * orhs_, i );
               row( osres_ , i ) = row( lhs_ * orhs_, i );
               row( refres_, i ) = row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( olhs_ * rhs_, i );
               row( odres_ , i ) = row( olhs_ * rhs_, i );
               row( sres_  , i ) = row( olhs_ * rhs_, i );
               row( osres_ , i ) = row( olhs_ * rhs_, i );
               row( refres_, i ) = row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( olhs_ * orhs_, i );
               row( odres_ , i ) = row( olhs_ * orhs_, i );
               row( sres_  , i ) = row( olhs_ * orhs_, i );
               row( osres_ , i ) = row( olhs_ * orhs_, i );
               row( refres_, i ) = row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Row-wise multiplication with evaluated matrices
      {
         test_  = "Row-wise multiplication with evaluated matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( eval( lhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) = row( eval( lhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) = row( eval( lhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) = row( eval( lhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) = row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( eval( lhs_ ) * eval( orhs_ ), i );
               row( odres_ , i ) = row( eval( lhs_ ) * eval( orhs_ ), i );
               row( sres_  , i ) = row( eval( lhs_ ) * eval( orhs_ ), i );
               row( osres_ , i ) = row( eval( lhs_ ) * eval( orhs_ ), i );
               row( refres_, i ) = row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( eval( olhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) = row( eval( olhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) = row( eval( olhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) = row( eval( olhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) = row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) = row( eval( olhs_ ) * eval( orhs_ ), i );
               row( odres_ , i ) = row( eval( olhs_ ) * eval( orhs_ ), i );
               row( sres_  , i ) = row( eval( olhs_ ) * eval( orhs_ ), i );
               row( osres_ , i ) = row( eval( olhs_ ) * eval( orhs_ ), i );
               row( refres_, i ) = row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Row-wise multiplication with addition assignment
      //=====================================================================================

      // Row-wise multiplication with addition assignment with the given matrices
      {
         test_  = "Row-wise multiplication with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( lhs_ * rhs_, i );
               row( odres_ , i ) += row( lhs_ * rhs_, i );
               row( sres_  , i ) += row( lhs_ * rhs_, i );
               row( osres_ , i ) += row( lhs_ * rhs_, i );
               row( refres_, i ) += row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( lhs_ * orhs_, i );
               row( odres_ , i ) += row( lhs_ * orhs_, i );
               row( sres_  , i ) += row( lhs_ * orhs_, i );
               row( osres_ , i ) += row( lhs_ * orhs_, i );
               row( refres_, i ) += row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( olhs_ * rhs_, i );
               row( odres_ , i ) += row( olhs_ * rhs_, i );
               row( sres_  , i ) += row( olhs_ * rhs_, i );
               row( osres_ , i ) += row( olhs_ * rhs_, i );
               row( refres_, i ) += row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( olhs_ * orhs_, i );
               row( odres_ , i ) += row( olhs_ * orhs_, i );
               row( sres_  , i ) += row( olhs_ * orhs_, i );
               row( osres_ , i ) += row( olhs_ * orhs_, i );
               row( refres_, i ) += row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Row-wise multiplication with addition assignment with evaluated matrices
      {
         test_  = "Row-wise multiplication with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( eval( lhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) += row( eval( lhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) += row( eval( lhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) += row( eval( lhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) += row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( eval( lhs_ ) * eval( orhs_ ), i );
               row( odres_ , i ) += row( eval( lhs_ ) * eval( orhs_ ), i );
               row( sres_  , i ) += row( eval( lhs_ ) * eval( orhs_ ), i );
               row( osres_ , i ) += row( eval( lhs_ ) * eval( orhs_ ), i );
               row( refres_, i ) += row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( eval( olhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) += row( eval( olhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) += row( eval( olhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) += row( eval( olhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) += row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) += row( eval( olhs_ ) * eval( orhs_ ), i );
               row( odres_ , i ) += row( eval( olhs_ ) * eval( orhs_ ), i );
               row( sres_  , i ) += row( eval( olhs_ ) * eval( orhs_ ), i );
               row( osres_ , i ) += row( eval( olhs_ ) * eval( orhs_ ), i );
               row( refres_, i ) += row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Row-wise multiplication with subtraction assignment
      //=====================================================================================

      // Row-wise multiplication with subtraction assignment with the given matrices
      {
         test_  = "Row-wise multiplication with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( lhs_ * rhs_, i );
               row( odres_ , i ) -= row( lhs_ * rhs_, i );
               row( sres_  , i ) -= row( lhs_ * rhs_, i );
               row( osres_ , i ) -= row( lhs_ * rhs_, i );
               row( refres_, i ) -= row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( lhs_ * orhs_, i );
               row( odres_ , i ) -= row( lhs_ * orhs_, i );
               row( sres_  , i ) -= row( lhs_ * orhs_, i );
               row( osres_ , i ) -= row( lhs_ * orhs_, i );
               row( refres_, i ) -= row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( olhs_ * rhs_, i );
               row( odres_ , i ) -= row( olhs_ * rhs_, i );
               row( sres_  , i ) -= row( olhs_ * rhs_, i );
               row( osres_ , i ) -= row( olhs_ * rhs_, i );
               row( refres_, i ) -= row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( olhs_ * orhs_, i );
               row( odres_ , i ) -= row( olhs_ * orhs_, i );
               row( sres_  , i ) -= row( olhs_ * orhs_, i );
               row( osres_ , i ) -= row( olhs_ * orhs_, i );
               row( refres_, i ) -= row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Row-wise multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Row-wise multiplication with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) -= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) -= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) -= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) -= row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( eval( lhs_ ) * eval( orhs_ ), i );
               row( odres_ , i ) -= row( eval( lhs_ ) * eval( orhs_ ), i );
               row( sres_  , i ) -= row( eval( lhs_ ) * eval( orhs_ ), i );
               row( osres_ , i ) -= row( eval( lhs_ ) * eval( orhs_ ), i );
               row( refres_, i ) -= row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( eval( olhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) -= row( eval( olhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) -= row( eval( olhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) -= row( eval( olhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) -= row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) -= row( eval( olhs_ ) * eval( orhs_ ), i );
               row( odres_ , i ) -= row( eval( olhs_ ) * eval( orhs_ ), i );
               row( sres_  , i ) -= row( eval( olhs_ ) * eval( orhs_ ), i );
               row( osres_ , i ) -= row( eval( olhs_ ) * eval( orhs_ ), i );
               row( refres_, i ) -= row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Row-wise multiplication with multiplication assignment
      //=====================================================================================

      // Row-wise multiplication with multiplication assignment with the given matrices
      {
         test_  = "Row-wise multiplication with multiplication assignment with the given matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( lhs_ * rhs_, i );
               row( odres_ , i ) *= row( lhs_ * rhs_, i );
               row( sres_  , i ) *= row( lhs_ * rhs_, i );
               row( osres_ , i ) *= row( lhs_ * rhs_, i );
               row( refres_, i ) *= row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( lhs_ * orhs_, i );
               row( odres_ , i ) *= row( lhs_ * orhs_, i );
               row( sres_  , i ) *= row( lhs_ * orhs_, i );
               row( osres_ , i ) *= row( lhs_ * orhs_, i );
               row( refres_, i ) *= row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( olhs_ * rhs_, i );
               row( odres_ , i ) *= row( olhs_ * rhs_, i );
               row( sres_  , i ) *= row( olhs_ * rhs_, i );
               row( osres_ , i ) *= row( olhs_ * rhs_, i );
               row( refres_, i ) *= row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( olhs_ * orhs_, i );
               row( odres_ , i ) *= row( olhs_ * orhs_, i );
               row( sres_  , i ) *= row( olhs_ * orhs_, i );
               row( osres_ , i ) *= row( olhs_ * orhs_, i );
               row( refres_, i ) *= row( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Row-wise multiplication with multiplication assignment with evaluated matrices
      {
         test_  = "Row-wise multiplication with multiplication assignment with evaluated matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) *= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) *= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) *= row( eval( lhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) *= row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( eval( lhs_ ) * eval( orhs_ ), i );
               row( odres_ , i ) *= row( eval( lhs_ ) * eval( orhs_ ), i );
               row( sres_  , i ) *= row( eval( lhs_ ) * eval( orhs_ ), i );
               row( osres_ , i ) *= row( eval( lhs_ ) * eval( orhs_ ), i );
               row( refres_, i ) *= row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( eval( olhs_ ) * eval( rhs_ ), i );
               row( odres_ , i ) *= row( eval( olhs_ ) * eval( rhs_ ), i );
               row( sres_  , i ) *= row( eval( olhs_ ) * eval( rhs_ ), i );
               row( osres_ , i ) *= row( eval( olhs_ ) * eval( rhs_ ), i );
               row( refres_, i ) *= row( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=0UL; i<lhs_.rows(); ++i ) {
               row( dres_  , i ) *= row( eval( olhs_ ) * eval( orhs_ ), i );
               row( odres_ , i ) *= row( eval( olhs_ ) * eval( orhs_ ), i );
               row( sres_  , i ) *= row( eval( olhs_ ) * eval( orhs_ ), i );
               row( osres_ , i ) *= row( eval( olhs_ ) * eval( orhs_ ), i );
               row( refres_, i ) *= row( eval( reflhs_ ) * eval( refrhs_ ), i );
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
            auto r = row( lhs_ * rhs_, lhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = row( lhs_ * orhs_, lhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = row( olhs_ * rhs_, olhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = row( olhs_ * orhs_, olhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
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
/*!\brief Testing the rows-wise sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the rows-wise matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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
      // Rows-wise multiplication
      //=====================================================================================

      // Rows-wise multiplication with the given matrices
      {
         test_  = "Rows-wise multiplication with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( lhs_ * orhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( lhs_ * orhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( lhs_ * orhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( lhs_ * orhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( reflhs_ * refrhs_, &indices[index], n );
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
               rows( dres_  , &indices[index], n ) = rows( olhs_ * rhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( olhs_ * rhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( olhs_ * rhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( olhs_ * rhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( reflhs_ * refrhs_, &indices[index], n );
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
               rows( dres_  , &indices[index], n ) = rows( olhs_ * orhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( olhs_ * orhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( olhs_ * orhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( olhs_ * orhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Rows-wise multiplication with evaluated matrices
      {
         test_  = "Rows-wise multiplication with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) = rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) = rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Rows-wise multiplication with addition assignment
      //=====================================================================================

      // Rows-wise multiplication with addition assignment with the given matrices
      {
         test_  = "Rows-wise multiplication with addition assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( lhs_ * orhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( lhs_ * orhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( lhs_ * orhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( lhs_ * orhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( reflhs_ * refrhs_, &indices[index], n );
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
               rows( dres_  , &indices[index], n ) += rows( olhs_ * rhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( olhs_ * rhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( olhs_ * rhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( olhs_ * rhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( reflhs_ * refrhs_, &indices[index], n );
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
               rows( dres_  , &indices[index], n ) += rows( olhs_ * orhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( olhs_ * orhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( olhs_ * orhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( olhs_ * orhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Rows-wise multiplication with addition assignment with evaluated matrices
      {
         test_  = "Rows-wise multiplication with addition assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) += rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) += rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Rows-wise multiplication with subtraction assignment
      //=====================================================================================

      // Rows-wise multiplication with subtraction assignment with the given matrices
      {
         test_  = "Rows-wise multiplication with subtraction assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( lhs_ * orhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( lhs_ * orhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( lhs_ * orhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( lhs_ * orhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( reflhs_ * refrhs_, &indices[index], n );
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
               rows( dres_  , &indices[index], n ) -= rows( olhs_ * rhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( olhs_ * rhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( olhs_ * rhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( olhs_ * rhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( reflhs_ * refrhs_, &indices[index], n );
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
               rows( dres_  , &indices[index], n ) -= rows( olhs_ * orhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( olhs_ * orhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( olhs_ * orhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( olhs_ * orhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Rows-wise multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Rows-wise multiplication with subtraction assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) -= rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) -= rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Rows-wise multiplication with Schur product assignment
      //=====================================================================================

      // Rows-wise multiplication with Schur product assignment with the given matrices
      {
         test_  = "Rows-wise multiplication with Schur product assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( lhs_ * orhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( lhs_ * orhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( lhs_ * orhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( lhs_ * orhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( reflhs_ * refrhs_, &indices[index], n );
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
               rows( dres_  , &indices[index], n ) %= rows( olhs_ * rhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( olhs_ * rhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( olhs_ * rhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( olhs_ * rhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( reflhs_ * refrhs_, &indices[index], n );
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
               rows( dres_  , &indices[index], n ) %= rows( olhs_ * orhs_, &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( olhs_ * orhs_, &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( olhs_ * orhs_, &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( olhs_ * orhs_, &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Rows-wise multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Rows-wise multiplication with Schur product assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) %= rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) %= rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
            auto r = rows( lhs_ * rhs_, { lhs_.rows() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( lhs_ * orhs_, { lhs_.rows() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( olhs_ * rhs_, { olhs_.rows() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( olhs_ * orhs_, { olhs_.rows() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
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
            auto r = rows( lhs_ * rhs_, [index=lhs_.rows()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( lhs_ * orhs_, [index=lhs_.rows()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( olhs_ * rhs_, [index=olhs_.rows()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( olhs_ * orhs_, [index=olhs_.rows()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
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
/*!\brief Skipping the rows-wise sparse matrix/sparse matrix multiplication.
//
// \return void
//
// This function is called in case the rows-wise matrix/matrix multiplication operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testRowsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the column-wise sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the column-wise matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testColumnOperation()
{
#if BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION
   if( BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION > 1 )
   {
      if( lhs_.columns() == 0UL )
         return;


      //=====================================================================================
      // Column-wise multiplication
      //=====================================================================================

      // Column-wise multiplication with the given matrices
      {
         test_  = "Column-wise multiplication with the given matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) = column( lhs_ * rhs_, j );
               column( odres_ , j ) = column( lhs_ * rhs_, j );
               column( sres_  , j ) = column( lhs_ * rhs_, j );
               column( osres_ , j ) = column( lhs_ * rhs_, j );
               column( refres_, j ) = column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         return;

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) = column( lhs_ * orhs_, j );
               column( odres_ , j ) = column( lhs_ * orhs_, j );
               column( sres_  , j ) = column( lhs_ * orhs_, j );
               column( osres_ , j ) = column( lhs_ * orhs_, j );
               column( refres_, j ) = column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) = column( olhs_ * rhs_, j );
               column( odres_ , j ) = column( olhs_ * rhs_, j );
               column( sres_  , j ) = column( olhs_ * rhs_, j );
               column( osres_ , j ) = column( olhs_ * rhs_, j );
               column( refres_, j ) = column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) = column( olhs_ * orhs_, j );
               column( odres_ , j ) = column( olhs_ * orhs_, j );
               column( sres_  , j ) = column( olhs_ * orhs_, j );
               column( osres_ , j ) = column( olhs_ * orhs_, j );
               column( refres_, j ) = column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Column-wise multiplication with evaluated matrices
      {
         test_  = "Column-wise multiplication with evaluated matrices";
         error_ = "Failed multiplication operation";

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) = column( eval( lhs_ ) * eval( rhs_ ), j );
               column( odres_ , j ) = column( eval( lhs_ ) * eval( rhs_ ), j );
               column( sres_  , j ) = column( eval( lhs_ ) * eval( rhs_ ), j );
               column( osres_ , j ) = column( eval( lhs_ ) * eval( rhs_ ), j );
               column( refres_, j ) = column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) = column( eval( lhs_ ) * eval( orhs_ ), j );
               column( odres_ , j ) = column( eval( lhs_ ) * eval( orhs_ ), j );
               column( sres_  , j ) = column( eval( lhs_ ) * eval( orhs_ ), j );
               column( osres_ , j ) = column( eval( lhs_ ) * eval( orhs_ ), j );
               column( refres_, j ) = column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) = column( eval( olhs_ ) * eval( rhs_ ), j );
               column( odres_ , j ) = column( eval( olhs_ ) * eval( rhs_ ), j );
               column( sres_  , j ) = column( eval( olhs_ ) * eval( rhs_ ), j );
               column( osres_ , j ) = column( eval( olhs_ ) * eval( rhs_ ), j );
               column( refres_, j ) = column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) = column( eval( olhs_ ) * eval( orhs_ ), j );
               column( odres_ , j ) = column( eval( olhs_ ) * eval( orhs_ ), j );
               column( sres_  , j ) = column( eval( olhs_ ) * eval( orhs_ ), j );
               column( osres_ , j ) = column( eval( olhs_ ) * eval( orhs_ ), j );
               column( refres_, j ) = column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Column-wise multiplication with addition assignment
      //=====================================================================================

      // Column-wise multiplication with addition assignment with the given matrices
      {
         test_  = "Column-wise multiplication with addition assignment with the given matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) += column( lhs_ * rhs_, j );
               column( odres_ , j ) += column( lhs_ * rhs_, j );
               column( sres_  , j ) += column( lhs_ * rhs_, j );
               column( osres_ , j ) += column( lhs_ * rhs_, j );
               column( refres_, j ) += column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) += column( lhs_ * orhs_, j );
               column( odres_ , j ) += column( lhs_ * orhs_, j );
               column( sres_  , j ) += column( lhs_ * orhs_, j );
               column( osres_ , j ) += column( lhs_ * orhs_, j );
               column( refres_, j ) += column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) += column( olhs_ * rhs_, j );
               column( odres_ , j ) += column( olhs_ * rhs_, j );
               column( sres_  , j ) += column( olhs_ * rhs_, j );
               column( osres_ , j ) += column( olhs_ * rhs_, j );
               column( refres_, j ) += column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) += column( olhs_ * orhs_, j );
               column( odres_ , j ) += column( olhs_ * orhs_, j );
               column( sres_  , j ) += column( olhs_ * orhs_, j );
               column( osres_ , j ) += column( olhs_ * orhs_, j );
               column( refres_, j ) += column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Column-wise multiplication with addition assignment with evaluated matrices
      {
         test_  = "Column-wise multiplication with addition assignment with evaluated matrices";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) += column( eval( lhs_ ) * eval( rhs_ ), j );
               column( odres_ , j ) += column( eval( lhs_ ) * eval( rhs_ ), j );
               column( sres_  , j ) += column( eval( lhs_ ) * eval( rhs_ ), j );
               column( osres_ , j ) += column( eval( lhs_ ) * eval( rhs_ ), j );
               column( refres_, j ) += column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) += column( eval( lhs_ ) * eval( orhs_ ), j );
               column( odres_ , j ) += column( eval( lhs_ ) * eval( orhs_ ), j );
               column( sres_  , j ) += column( eval( lhs_ ) * eval( orhs_ ), j );
               column( osres_ , j ) += column( eval( lhs_ ) * eval( orhs_ ), j );
               column( refres_, j ) += column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) += column( eval( olhs_ ) * eval( rhs_ ), j );
               column( odres_ , j ) += column( eval( olhs_ ) * eval( rhs_ ), j );
               column( sres_  , j ) += column( eval( olhs_ ) * eval( rhs_ ), j );
               column( osres_ , j ) += column( eval( olhs_ ) * eval( rhs_ ), j );
               column( refres_, j ) += column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) += column( eval( olhs_ ) * eval( orhs_ ), j );
               column( odres_ , j ) += column( eval( olhs_ ) * eval( orhs_ ), j );
               column( sres_  , j ) += column( eval( olhs_ ) * eval( orhs_ ), j );
               column( osres_ , j ) += column( eval( olhs_ ) * eval( orhs_ ), j );
               column( refres_, j ) += column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Column-wise multiplication with subtraction assignment
      //=====================================================================================

      // Column-wise multiplication with subtraction assignment with the given matrices
      {
         test_  = "Column-wise multiplication with subtraction assignment with the given matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( lhs_ * rhs_, j );
               column( odres_ , j ) -= column( lhs_ * rhs_, j );
               column( sres_  , j ) -= column( lhs_ * rhs_, j );
               column( osres_ , j ) -= column( lhs_ * rhs_, j );
               column( refres_, j ) -= column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( lhs_ * orhs_, j );
               column( odres_ , j ) -= column( lhs_ * orhs_, j );
               column( sres_  , j ) -= column( lhs_ * orhs_, j );
               column( osres_ , j ) -= column( lhs_ * orhs_, j );
               column( refres_, j ) -= column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( olhs_ * rhs_, j );
               column( odres_ , j ) -= column( olhs_ * rhs_, j );
               column( sres_  , j ) -= column( olhs_ * rhs_, j );
               column( osres_ , j ) -= column( olhs_ * rhs_, j );
               column( refres_, j ) -= column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( olhs_ * orhs_, j );
               column( odres_ , j ) -= column( olhs_ * orhs_, j );
               column( sres_  , j ) -= column( olhs_ * orhs_, j );
               column( osres_ , j ) -= column( olhs_ * orhs_, j );
               column( refres_, j ) -= column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Column-wise multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Column-wise multiplication with subtraction assignment with evaluated matrices";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( eval( lhs_ ) * eval( rhs_ ), j );
               column( odres_ , j ) -= column( eval( lhs_ ) * eval( rhs_ ), j );
               column( sres_  , j ) -= column( eval( lhs_ ) * eval( rhs_ ), j );
               column( osres_ , j ) -= column( eval( lhs_ ) * eval( rhs_ ), j );
               column( refres_, j ) -= column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( eval( lhs_ ) * eval( orhs_ ), j );
               column( odres_ , j ) -= column( eval( lhs_ ) * eval( orhs_ ), j );
               column( sres_  , j ) -= column( eval( lhs_ ) * eval( orhs_ ), j );
               column( osres_ , j ) -= column( eval( lhs_ ) * eval( orhs_ ), j );
               column( refres_, j ) -= column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( eval( olhs_ ) * eval( rhs_ ), j );
               column( odres_ , j ) -= column( eval( olhs_ ) * eval( rhs_ ), j );
               column( sres_  , j ) -= column( eval( olhs_ ) * eval( rhs_ ), j );
               column( osres_ , j ) -= column( eval( olhs_ ) * eval( rhs_ ), j );
               column( refres_, j ) -= column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) -= column( eval( olhs_ ) * eval( orhs_ ), j );
               column( odres_ , j ) -= column( eval( olhs_ ) * eval( orhs_ ), j );
               column( sres_  , j ) -= column( eval( olhs_ ) * eval( orhs_ ), j );
               column( osres_ , j ) -= column( eval( olhs_ ) * eval( orhs_ ), j );
               column( refres_, j ) -= column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Column-wise multiplication with multiplication assignment
      //=====================================================================================

      // Column-wise multiplication with multiplication assignment with the given matrices
      {
         test_  = "Column-wise multiplication with multiplication assignment with the given matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( lhs_ * rhs_, j );
               column( odres_ , j ) *= column( lhs_ * rhs_, j );
               column( sres_  , j ) *= column( lhs_ * rhs_, j );
               column( osres_ , j ) *= column( lhs_ * rhs_, j );
               column( refres_, j ) *= column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( lhs_ * orhs_, j );
               column( odres_ , j ) *= column( lhs_ * orhs_, j );
               column( sres_  , j ) *= column( lhs_ * orhs_, j );
               column( osres_ , j ) *= column( lhs_ * orhs_, j );
               column( refres_, j ) *= column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( olhs_ * rhs_, j );
               column( odres_ , j ) *= column( olhs_ * rhs_, j );
               column( sres_  , j ) *= column( olhs_ * rhs_, j );
               column( osres_ , j ) *= column( olhs_ * rhs_, j );
               column( refres_, j ) *= column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( olhs_ * orhs_, j );
               column( odres_ , j ) *= column( olhs_ * orhs_, j );
               column( sres_  , j ) *= column( olhs_ * orhs_, j );
               column( osres_ , j ) *= column( olhs_ * orhs_, j );
               column( refres_, j ) *= column( reflhs_ * refrhs_, j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Column-wise multiplication with multiplication assignment with evaluated matrices
      {
         test_  = "Column-wise multiplication with multiplication assignment with evaluated matrices";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( eval( lhs_ ) * eval( rhs_ ), j );
               column( odres_ , j ) *= column( eval( lhs_ ) * eval( rhs_ ), j );
               column( sres_  , j ) *= column( eval( lhs_ ) * eval( rhs_ ), j );
               column( osres_ , j ) *= column( eval( lhs_ ) * eval( rhs_ ), j );
               column( refres_, j ) *= column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( eval( lhs_ ) * eval( orhs_ ), j );
               column( odres_ , j ) *= column( eval( lhs_ ) * eval( orhs_ ), j );
               column( sres_  , j ) *= column( eval( lhs_ ) * eval( orhs_ ), j );
               column( osres_ , j ) *= column( eval( lhs_ ) * eval( orhs_ ), j );
               column( refres_, j ) *= column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( eval( olhs_ ) * eval( rhs_ ), j );
               column( odres_ , j ) *= column( eval( olhs_ ) * eval( rhs_ ), j );
               column( sres_  , j ) *= column( eval( olhs_ ) * eval( rhs_ ), j );
               column( osres_ , j ) *= column( eval( olhs_ ) * eval( rhs_ ), j );
               column( refres_, j ) *= column( eval( reflhs_ ) * eval( refrhs_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t j=0UL; j<rhs_.columns(); ++j ) {
               column( dres_  , j ) *= column( eval( olhs_ ) * eval( orhs_ ), j );
               column( odres_ , j ) *= column( eval( olhs_ ) * eval( orhs_ ), j );
               column( sres_  , j ) *= column( eval( olhs_ ) * eval( orhs_ ), j );
               column( osres_ , j ) *= column( eval( olhs_ ) * eval( orhs_ ), j );
               column( refres_, j ) *= column( eval( reflhs_ ) * eval( refrhs_ ), j );
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
            auto c = column( lhs_ * rhs_, rhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = column( lhs_ * orhs_, orhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = column( olhs_ * rhs_, rhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = column( olhs_ * orhs_, orhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
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
/*!\brief Testing the columns-wise sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the columns-wise matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testColumnsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION > 1 )
   {
      if( rhs_.columns() == 0UL )
         return;


      std::vector<size_t> indices( rhs_.columns() );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );


      //=====================================================================================
      // Columns-wise multiplication
      //=====================================================================================

      // Columns-wise multiplication with the given matrices
      {
         test_  = "Columns-wise multiplication with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( lhs_ * orhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( lhs_ * orhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( lhs_ * orhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( lhs_ * orhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( reflhs_ * refrhs_, &indices[index], n );
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
               columns( dres_  , &indices[index], n ) = columns( olhs_ * rhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( olhs_ * rhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( olhs_ * rhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( olhs_ * rhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( reflhs_ * refrhs_, &indices[index], n );
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
               columns( dres_  , &indices[index], n ) = columns( olhs_ * orhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( olhs_ * orhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( olhs_ * orhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( olhs_ * orhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Columns-wise multiplication with evaluated matrices
      {
         test_  = "Columns-wise multiplication with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) = columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) = columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Columns-wise multiplication with addition assignment
      //=====================================================================================

      // Columns-wise multiplication with addition assignment with the given matrices
      {
         test_  = "Columns-wise multiplication with addition assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( lhs_ * orhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( lhs_ * orhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( lhs_ * orhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( lhs_ * orhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( reflhs_ * refrhs_, &indices[index], n );
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
               columns( dres_  , &indices[index], n ) += columns( olhs_ * rhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( olhs_ * rhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( olhs_ * rhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( olhs_ * rhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( reflhs_ * refrhs_, &indices[index], n );
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
               columns( dres_  , &indices[index], n ) += columns( olhs_ * orhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( olhs_ * orhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( olhs_ * orhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( olhs_ * orhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Columns-wise multiplication with addition assignment with evaluated matrices
      {
         test_  = "Columns-wise multiplication with addition assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) += columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) += columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Columns-wise multiplication with subtraction assignment
      //=====================================================================================

      // Columns-wise multiplication with subtraction assignment with the given matrices
      {
         test_  = "Columns-wise multiplication with subtraction assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( lhs_ * orhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( lhs_ * orhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( lhs_ * orhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( lhs_ * orhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( reflhs_ * refrhs_, &indices[index], n );
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
               columns( dres_  , &indices[index], n ) -= columns( olhs_ * rhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( olhs_ * rhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( olhs_ * rhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( olhs_ * rhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( reflhs_ * refrhs_, &indices[index], n );
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
               columns( dres_  , &indices[index], n ) -= columns( olhs_ * orhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( olhs_ * orhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( olhs_ * orhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( olhs_ * orhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Columns-wise multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Columns-wise multiplication with subtraction assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) -= columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) -= columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Columns-wise multiplication with Schur product assignment
      //=====================================================================================

      // Columns-wise multiplication with Schur product assignment with the given matrices
      {
         test_  = "Columns-wise multiplication with Schur product assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( lhs_ * orhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( lhs_ * orhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( lhs_ * orhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( lhs_ * orhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( reflhs_ * refrhs_, &indices[index], n );
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
               columns( dres_  , &indices[index], n ) %= columns( olhs_ * rhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( olhs_ * rhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( olhs_ * rhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( olhs_ * rhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( reflhs_ * refrhs_, &indices[index], n );
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
               columns( dres_  , &indices[index], n ) %= columns( olhs_ * orhs_, &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( olhs_ * orhs_, &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( olhs_ * orhs_, &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( olhs_ * orhs_, &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( reflhs_ * refrhs_, &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Columns-wise multiplication with Schur product assignment with evaluated matrices
      {
         test_  = "Columns-wise multiplication with Schur product assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( eval( lhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) %= columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( eval( olhs_ ) * eval( rhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) %= columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( eval( olhs_ ) * eval( orhs_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( eval( reflhs_ ) * eval( refrhs_ ), &indices[index], n );
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
            auto c = columns( lhs_ * rhs_, { rhs_.columns() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( lhs_ * orhs_, { orhs_.columns() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( olhs_ * rhs_, { rhs_.columns() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( olhs_ * orhs_, { orhs_.columns() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
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
            auto c = columns( lhs_ * rhs_, [index=rhs_.columns()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( lhs_ * orhs_, [index=orhs_.columns()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( olhs_ * rhs_, [index=rhs_.columns()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( olhs_ * orhs_, [index=orhs_.columns()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
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
/*!\brief Skipping the columns-wise sparse matrix/sparse matrix multiplication.
//
// \return void
//
// This function is called in case the columns-wise matrix/matrix multiplication operation is not
// available for the given matrix types \a MT1 and \a MT2.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testColumnsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the band-wise sparse matrix/sparse matrix multiplication.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the band-wise matrix multiplication with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the multiplication or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::testBandOperation()
{
#if BLAZETEST_MATHTEST_TEST_BAND_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BAND_OPERATION > 1 )
   {
      if( lhs_.rows() == 0UL || rhs_.columns() == 0UL )
         return;


      const ptrdiff_t ibegin( 1UL - lhs_.rows() );
      const ptrdiff_t iend  ( rhs_.columns() );


      //=====================================================================================
      // Band-wise multiplication
      //=====================================================================================

      // Band-wise multiplication with the given matrices
      {
         test_  = "Band-wise multiplication with the given matrices";
         error_ = "Failed multiplication operation";

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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( lhs_ * orhs_, i );
               band( odres_ , i ) = band( lhs_ * orhs_, i );
               band( sres_  , i ) = band( lhs_ * orhs_, i );
               band( osres_ , i ) = band( lhs_ * orhs_, i );
               band( refres_, i ) = band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( olhs_ * rhs_, i );
               band( odres_ , i ) = band( olhs_ * rhs_, i );
               band( sres_  , i ) = band( olhs_ * rhs_, i );
               band( osres_ , i ) = band( olhs_ * rhs_, i );
               band( refres_, i ) = band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( size_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( olhs_ * orhs_, i );
               band( odres_ , i ) = band( olhs_ * orhs_, i );
               band( sres_  , i ) = band( olhs_ * orhs_, i );
               band( osres_ , i ) = band( olhs_ * orhs_, i );
               band( refres_, i ) = band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Band-wise multiplication with evaluated matrices
      {
         test_  = "Band-wise multiplication with evaluated matrices";
         error_ = "Failed multiplication operation";

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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( eval( lhs_ ) * eval( orhs_ ), i );
               band( odres_ , i ) = band( eval( lhs_ ) * eval( orhs_ ), i );
               band( sres_  , i ) = band( eval( lhs_ ) * eval( orhs_ ), i );
               band( osres_ , i ) = band( eval( lhs_ ) * eval( orhs_ ), i );
               band( refres_, i ) = band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( eval( olhs_ ) * eval( rhs_ ), i );
               band( odres_ , i ) = band( eval( olhs_ ) * eval( rhs_ ), i );
               band( sres_  , i ) = band( eval( olhs_ ) * eval( rhs_ ), i );
               band( osres_ , i ) = band( eval( olhs_ ) * eval( rhs_ ), i );
               band( refres_, i ) = band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( eval( olhs_ ) * eval( orhs_ ), i );
               band( odres_ , i ) = band( eval( olhs_ ) * eval( orhs_ ), i );
               band( sres_  , i ) = band( eval( olhs_ ) * eval( orhs_ ), i );
               band( osres_ , i ) = band( eval( olhs_ ) * eval( orhs_ ), i );
               band( refres_, i ) = band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Band-wise multiplication with addition assignment
      //=====================================================================================

      // Band-wise multiplication with addition assignment with the given matrices
      {
         test_  = "Band-wise multiplication with addition assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( lhs_ * orhs_, i );
               band( odres_ , i ) += band( lhs_ * orhs_, i );
               band( sres_  , i ) += band( lhs_ * orhs_, i );
               band( osres_ , i ) += band( lhs_ * orhs_, i );
               band( refres_, i ) += band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( olhs_ * rhs_, i );
               band( odres_ , i ) += band( olhs_ * rhs_, i );
               band( sres_  , i ) += band( olhs_ * rhs_, i );
               band( osres_ , i ) += band( olhs_ * rhs_, i );
               band( refres_, i ) += band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( olhs_ * orhs_, i );
               band( odres_ , i ) += band( olhs_ * orhs_, i );
               band( sres_  , i ) += band( olhs_ * orhs_, i );
               band( osres_ , i ) += band( olhs_ * orhs_, i );
               band( refres_, i ) += band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Band-wise multiplication with addition assignment with evaluated matrices
      {
         test_  = "Band-wise multiplication with addition assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( eval( lhs_ ) * eval( orhs_ ), i );
               band( odres_ , i ) += band( eval( lhs_ ) * eval( orhs_ ), i );
               band( sres_  , i ) += band( eval( lhs_ ) * eval( orhs_ ), i );
               band( osres_ , i ) += band( eval( lhs_ ) * eval( orhs_ ), i );
               band( refres_, i ) += band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( eval( olhs_ ) * eval( rhs_ ), i );
               band( odres_ , i ) += band( eval( olhs_ ) * eval( rhs_ ), i );
               band( sres_  , i ) += band( eval( olhs_ ) * eval( rhs_ ), i );
               band( osres_ , i ) += band( eval( olhs_ ) * eval( rhs_ ), i );
               band( refres_, i ) += band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( eval( olhs_ ) * eval( orhs_ ), i );
               band( odres_ , i ) += band( eval( olhs_ ) * eval( orhs_ ), i );
               band( sres_  , i ) += band( eval( olhs_ ) * eval( orhs_ ), i );
               band( osres_ , i ) += band( eval( olhs_ ) * eval( orhs_ ), i );
               band( refres_, i ) += band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Band-wise multiplication with subtraction assignment
      //=====================================================================================

      // Band-wise multiplication with subtraction assignment with the given matrices
      {
         test_  = "Band-wise multiplication with subtraction assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( lhs_ * orhs_, i );
               band( odres_ , i ) -= band( lhs_ * orhs_, i );
               band( sres_  , i ) -= band( lhs_ * orhs_, i );
               band( osres_ , i ) -= band( lhs_ * orhs_, i );
               band( refres_, i ) -= band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( olhs_ * rhs_, i );
               band( odres_ , i ) -= band( olhs_ * rhs_, i );
               band( sres_  , i ) -= band( olhs_ * rhs_, i );
               band( osres_ , i ) -= band( olhs_ * rhs_, i );
               band( refres_, i ) -= band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( olhs_ * orhs_, i );
               band( odres_ , i ) -= band( olhs_ * orhs_, i );
               band( sres_  , i ) -= band( olhs_ * orhs_, i );
               band( osres_ , i ) -= band( olhs_ * orhs_, i );
               band( refres_, i ) -= band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Band-wise multiplication with subtraction assignment with evaluated matrices
      {
         test_  = "Band-wise multiplication with subtraction assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( eval( lhs_ ) * eval( orhs_ ), i );
               band( odres_ , i ) -= band( eval( lhs_ ) * eval( orhs_ ), i );
               band( sres_  , i ) -= band( eval( lhs_ ) * eval( orhs_ ), i );
               band( osres_ , i ) -= band( eval( lhs_ ) * eval( orhs_ ), i );
               band( refres_, i ) -= band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( eval( olhs_ ) * eval( rhs_ ), i );
               band( odres_ , i ) -= band( eval( olhs_ ) * eval( rhs_ ), i );
               band( sres_  , i ) -= band( eval( olhs_ ) * eval( rhs_ ), i );
               band( osres_ , i ) -= band( eval( olhs_ ) * eval( rhs_ ), i );
               band( refres_, i ) -= band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( eval( olhs_ ) * eval( orhs_ ), i );
               band( odres_ , i ) -= band( eval( olhs_ ) * eval( orhs_ ), i );
               band( sres_  , i ) -= band( eval( olhs_ ) * eval( orhs_ ), i );
               band( osres_ , i ) -= band( eval( olhs_ ) * eval( orhs_ ), i );
               band( refres_, i ) -= band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }


      //=====================================================================================
      // Band-wise multiplication with multiplication assignment
      //=====================================================================================

      // Band-wise multiplication with multiplication assignment with the given matrices
      {
         test_  = "Band-wise multiplication with multiplication assignment with the given matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( lhs_ * orhs_, i );
               band( odres_ , i ) *= band( lhs_ * orhs_, i );
               band( sres_  , i ) *= band( lhs_ * orhs_, i );
               band( osres_ , i ) *= band( lhs_ * orhs_, i );
               band( refres_, i ) *= band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( olhs_ * rhs_, i );
               band( odres_ , i ) *= band( olhs_ * rhs_, i );
               band( sres_  , i ) *= band( olhs_ * rhs_, i );
               band( osres_ , i ) *= band( olhs_ * rhs_, i );
               band( refres_, i ) *= band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( olhs_ * orhs_, i );
               band( odres_ , i ) *= band( olhs_ * orhs_, i );
               band( sres_  , i ) *= band( olhs_ * orhs_, i );
               band( osres_ , i ) *= band( olhs_ * orhs_, i );
               band( refres_, i ) *= band( reflhs_ * refrhs_, i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,OMT2>( ex );
         }

         checkResults<OMT1,OMT2>();
      }

      // Band-wise multiplication with multiplication assignment with evaluated matrices
      {
         test_  = "Band-wise multiplication with multiplication assignment with evaluated matrices";
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
            convertException<MT1,MT2>( ex );
         }

         checkResults<MT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( eval( lhs_ ) * eval( orhs_ ), i );
               band( odres_ , i ) *= band( eval( lhs_ ) * eval( orhs_ ), i );
               band( sres_  , i ) *= band( eval( lhs_ ) * eval( orhs_ ), i );
               band( osres_ , i ) *= band( eval( lhs_ ) * eval( orhs_ ), i );
               band( refres_, i ) *= band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT1,OMT2>( ex );
         }

         checkResults<MT1,OMT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( eval( olhs_ ) * eval( rhs_ ), i );
               band( odres_ , i ) *= band( eval( olhs_ ) * eval( rhs_ ), i );
               band( sres_  , i ) *= band( eval( olhs_ ) * eval( rhs_ ), i );
               band( osres_ , i ) *= band( eval( olhs_ ) * eval( rhs_ ), i );
               band( refres_, i ) *= band( eval( reflhs_ ) * eval( refrhs_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT1,MT2>( ex );
         }

         checkResults<OMT1,MT2>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( eval( olhs_ ) * eval( orhs_ ), i );
               band( odres_ , i ) *= band( eval( olhs_ ) * eval( orhs_ ), i );
               band( sres_  , i ) *= band( eval( olhs_ ) * eval( orhs_ ), i );
               band( osres_ , i ) *= band( eval( olhs_ ) * eval( orhs_ ), i );
               band( refres_, i ) *= band( eval( reflhs_ ) * eval( refrhs_ ), i );
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
            auto b = band( lhs_ * rhs_, -lhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( lhs_ * orhs_, -lhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( olhs_ * rhs_, -olhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( olhs_ * orhs_, -olhs_.rows() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
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
            auto b = band( lhs_ * rhs_, rhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( lhs_ * orhs_, orhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( MT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( OMT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( olhs_ * rhs_, rhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
                << "     " << typeid( MT2 ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( olhs_ * orhs_, orhs_.columns() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Left-hand side sparse matrix type:\n"
                << "     " << typeid( OMT1 ).name() << "\n"
                << "   Right-hand side sparse matrix type:\n"
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
/*!\brief Testing the customized sparse matrix/sparse matrix multiplication.
//
// \param op The custom operation to be tested.
// \param name The human-readable name of the operation.
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the matrix multiplication with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment in combination with a custom operation.
// In case any error resulting from the multiplication or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
template< typename OP >   // Type of the custom operation
void OperationTest<MT1,MT2>::testCustomOperation( OP op, const std::string& name )
{
   //=====================================================================================
   // Customized multiplication
   //=====================================================================================

   // Customized multiplication with the given matrices
   {
      test_  = "Customized multiplication with the given matrices (" + name + ")";
      error_ = "Failed multiplication operation";

      try {
         initResults();
         dres_   = op( lhs_ * rhs_ );
         odres_  = op( lhs_ * rhs_ );
         sres_   = op( lhs_ * rhs_ );
         osres_  = op( lhs_ * rhs_ );
         refres_ = op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   = op( lhs_ * orhs_ );
         odres_  = op( lhs_ * orhs_ );
         sres_   = op( lhs_ * orhs_ );
         osres_  = op( lhs_ * orhs_ );
         refres_ = op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   = op( olhs_ * rhs_ );
         odres_  = op( olhs_ * rhs_ );
         sres_   = op( olhs_ * rhs_ );
         osres_  = op( olhs_ * rhs_ );
         refres_ = op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   = op( olhs_ * orhs_ );
         odres_  = op( olhs_ * orhs_ );
         sres_   = op( olhs_ * orhs_ );
         osres_  = op( olhs_ * orhs_ );
         refres_ = op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }

   // Customized multiplication with evaluated matrices
   {
      test_  = "Customized multiplication with evaluated matrices (" + name + ")";
      error_ = "Failed multiplication operation";

      try {
         initResults();
         dres_   = op( eval( lhs_ ) * eval( rhs_ ) );
         odres_  = op( eval( lhs_ ) * eval( rhs_ ) );
         sres_   = op( eval( lhs_ ) * eval( rhs_ ) );
         osres_  = op( eval( lhs_ ) * eval( rhs_ ) );
         refres_ = op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   = op( eval( lhs_ ) * eval( orhs_ ) );
         odres_  = op( eval( lhs_ ) * eval( orhs_ ) );
         sres_   = op( eval( lhs_ ) * eval( orhs_ ) );
         osres_  = op( eval( lhs_ ) * eval( orhs_ ) );
         refres_ = op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   = op( eval( olhs_ ) * eval( rhs_ ) );
         odres_  = op( eval( olhs_ ) * eval( rhs_ ) );
         sres_   = op( eval( olhs_ ) * eval( rhs_ ) );
         osres_  = op( eval( olhs_ ) * eval( rhs_ ) );
         refres_ = op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   = op( eval( olhs_ ) * eval( orhs_ ) );
         odres_  = op( eval( olhs_ ) * eval( orhs_ ) );
         sres_   = op( eval( olhs_ ) * eval( orhs_ ) );
         osres_  = op( eval( olhs_ ) * eval( orhs_ ) );
         refres_ = op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }


   //=====================================================================================
   // Customized multiplication with addition assignment
   //=====================================================================================

   // Customized multiplication with addition assignment with the given matrices
   {
      test_  = "Customized multiplication with addition assignment with the given matrices (" + name + ")";
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
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   += op( lhs_ * orhs_ );
         odres_  += op( lhs_ * orhs_ );
         sres_   += op( lhs_ * orhs_ );
         osres_  += op( lhs_ * orhs_ );
         refres_ += op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   += op( olhs_ * rhs_ );
         odres_  += op( olhs_ * rhs_ );
         sres_   += op( olhs_ * rhs_ );
         osres_  += op( olhs_ * rhs_ );
         refres_ += op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   += op( olhs_ * orhs_ );
         odres_  += op( olhs_ * orhs_ );
         sres_   += op( olhs_ * orhs_ );
         osres_  += op( olhs_ * orhs_ );
         refres_ += op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }

   // Customized multiplication with addition assignment with evaluated matrices
   {
      test_  = "Customized multiplication with addition assignment with evaluated matrices (" + name + ")";
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
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   += op( eval( lhs_ ) * eval( orhs_ ) );
         odres_  += op( eval( lhs_ ) * eval( orhs_ ) );
         sres_   += op( eval( lhs_ ) * eval( orhs_ ) );
         osres_  += op( eval( lhs_ ) * eval( orhs_ ) );
         refres_ += op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   += op( eval( olhs_ ) * eval( rhs_ ) );
         odres_  += op( eval( olhs_ ) * eval( rhs_ ) );
         sres_   += op( eval( olhs_ ) * eval( rhs_ ) );
         osres_  += op( eval( olhs_ ) * eval( rhs_ ) );
         refres_ += op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   += op( eval( olhs_ ) * eval( orhs_ ) );
         odres_  += op( eval( olhs_ ) * eval( orhs_ ) );
         sres_   += op( eval( olhs_ ) * eval( orhs_ ) );
         osres_  += op( eval( olhs_ ) * eval( orhs_ ) );
         refres_ += op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }


   //=====================================================================================
   // Customized multiplication with subtraction assignment
   //=====================================================================================

   // Customized multiplication with subtraction assignment with the given matrices
   {
      test_  = "Customized multiplication with subtraction assignment with the given matrices (" + name + ")";
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
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   -= op( lhs_ * orhs_ );
         odres_  -= op( lhs_ * orhs_ );
         sres_   -= op( lhs_ * orhs_ );
         osres_  -= op( lhs_ * orhs_ );
         refres_ -= op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   -= op( olhs_ * rhs_ );
         odres_  -= op( olhs_ * rhs_ );
         sres_   -= op( olhs_ * rhs_ );
         osres_  -= op( olhs_ * rhs_ );
         refres_ -= op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   -= op( olhs_ * orhs_ );
         odres_  -= op( olhs_ * orhs_ );
         sres_   -= op( olhs_ * orhs_ );
         osres_  -= op( olhs_ * orhs_ );
         refres_ -= op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }

   // Customized multiplication with subtraction assignment with evaluated matrices
   {
      test_  = "Customized multiplication with subtraction assignment with evaluated matrices (" + name + ")";
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
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   -= op( eval( lhs_ ) * eval( orhs_ ) );
         odres_  -= op( eval( lhs_ ) * eval( orhs_ ) );
         sres_   -= op( eval( lhs_ ) * eval( orhs_ ) );
         osres_  -= op( eval( lhs_ ) * eval( orhs_ ) );
         refres_ -= op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   -= op( eval( olhs_ ) * eval( rhs_ ) );
         odres_  -= op( eval( olhs_ ) * eval( rhs_ ) );
         sres_   -= op( eval( olhs_ ) * eval( rhs_ ) );
         osres_  -= op( eval( olhs_ ) * eval( rhs_ ) );
         refres_ -= op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   -= op( eval( olhs_ ) * eval( orhs_ ) );
         odres_  -= op( eval( olhs_ ) * eval( orhs_ ) );
         sres_   -= op( eval( olhs_ ) * eval( orhs_ ) );
         osres_  -= op( eval( olhs_ ) * eval( orhs_ ) );
         refres_ -= op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }


   //=====================================================================================
   // Customized multiplication with Schur product assignment
   //=====================================================================================

   // Customized multiplication with Schur product assignment with the given matrices
   {
      test_  = "Customized multiplication with Schur product assignment with the given matrices (" + name + ")";
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
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   %= op( lhs_ * orhs_ );
         odres_  %= op( lhs_ * orhs_ );
         sres_   %= op( lhs_ * orhs_ );
         osres_  %= op( lhs_ * orhs_ );
         refres_ %= op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   %= op( olhs_ * rhs_ );
         odres_  %= op( olhs_ * rhs_ );
         sres_   %= op( olhs_ * rhs_ );
         osres_  %= op( olhs_ * rhs_ );
         refres_ %= op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   %= op( olhs_ * orhs_ );
         odres_  %= op( olhs_ * orhs_ );
         sres_   %= op( olhs_ * orhs_ );
         osres_  %= op( olhs_ * orhs_ );
         refres_ %= op( reflhs_ * refrhs_ );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,OMT2>( ex );
      }

      checkResults<OMT1,OMT2>();
   }

   // Customized multiplication with Schur product assignment with evaluated matrices
   {
      test_  = "Customized multiplication with Schur product assignment with evaluated matrices (" + name + ")";
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
         convertException<MT1,MT2>( ex );
      }

      checkResults<MT1,MT2>();

      try {
         initResults();
         dres_   %= op( eval( lhs_ ) * eval( orhs_ ) );
         odres_  %= op( eval( lhs_ ) * eval( orhs_ ) );
         sres_   %= op( eval( lhs_ ) * eval( orhs_ ) );
         osres_  %= op( eval( lhs_ ) * eval( orhs_ ) );
         refres_ %= op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT1,OMT2>( ex );
      }

      checkResults<MT1,OMT2>();

      try {
         initResults();
         dres_   %= op( eval( olhs_ ) * eval( rhs_ ) );
         odres_  %= op( eval( olhs_ ) * eval( rhs_ ) );
         sres_   %= op( eval( olhs_ ) * eval( rhs_ ) );
         osres_  %= op( eval( olhs_ ) * eval( rhs_ ) );
         refres_ %= op( eval( reflhs_ ) * eval( refrhs_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT1,MT2>( ex );
      }

      checkResults<OMT1,MT2>();

      try {
         initResults();
         dres_   %= op( eval( olhs_ ) * eval( orhs_ ) );
         odres_  %= op( eval( olhs_ ) * eval( orhs_ ) );
         sres_   %= op( eval( olhs_ ) * eval( orhs_ ) );
         osres_  %= op( eval( olhs_ ) * eval( orhs_ ) );
         refres_ %= op( eval( reflhs_ ) * eval( refrhs_ ) );
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
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
          << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
          << "     " << typeid( LT ).name() << "\n"
          << "   Right-hand side " << ( IsRowMajorMatrix<RT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::initResults()
{
   const blaze::UnderlyingBuiltin_t<SRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<SRE> max( randmax );

   resize( sres_, rows( lhs_ ), columns( rhs_ ) );
   randomize( sres_, min, max );

   dres_   = sres_;
   odres_  = sres_;
   osres_  = sres_;
   refres_ = sres_;
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
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void OperationTest<MT1,MT2>::initTransposeResults()
{
   const blaze::UnderlyingBuiltin_t<TSRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<TSRE> max( randmax );

   resize( tsres_, columns( rhs_ ), rows( lhs_ ) );
   randomize( tsres_, min, max );

   tdres_  = tsres_;
   todres_ = tsres_;
   tosres_ = tsres_;
   refres_ = tsres_;
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
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
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
       << "   Left-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
       << "     " << typeid( LT ).name() << "\n"
       << "   Right-hand side " << ( IsRowMajorMatrix<LT>::value ? ( "row-major" ) : ( "column-major" ) ) << " sparse matrix type:\n"
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
/*!\brief Testing the matrix multiplication between two specific matrix types.
//
// \param creator1 The creator for the left-hand side matrix.
// \param creator2 The creator for the right-hand side matrix.
// \return void
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
void runTest( const Creator<MT1>& creator1, const Creator<MT2>& creator2 )
{
#if BLAZETEST_MATHTEST_TEST_MULTIPLICATION
   if( BLAZETEST_MATHTEST_TEST_MULTIPLICATION > 1 )
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
/*!\brief Macro for the definition of a sparse matrix/sparse matrix multiplication test case.
*/
#define DEFINE_SMATSMATMULT_OPERATION_TEST( MT1, MT2 ) \
   extern template class blazetest::mathtest::operations::smatsmatmult::OperationTest<MT1,MT2>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a sparse matrix/sparse matrix multiplication test case.
*/
#define RUN_SMATSMATMULT_OPERATION_TEST( C1, C2 ) \
   blazetest::mathtest::operations::smatsmatmult::runTest( C1, C2 )
/*! \endcond */
//*************************************************************************************************

} // namespace smatsmatmult

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
