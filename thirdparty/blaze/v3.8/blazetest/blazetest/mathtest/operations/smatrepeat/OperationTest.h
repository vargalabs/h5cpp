//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/smatrepeat/OperationTest.h
//  \brief Header file for the sparse matrix repeat operation test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_SMATREPEAT_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_SMATREPEAT_OPERATIONTEST_H_


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
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Functors.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/traits/RepeatTrait.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingScalar.h>
#include <blaze/math/Views.h>
#include <blaze/util/constraints/DerivedFrom.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/Nor.h>
#include <blaze/util/mpl/Not.h>
#include <blaze/util/Random.h>
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

namespace smatrepeat {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the sparse matrix repeat operation test.
//
// This class template represents one particular test of a repeat operation on a matrix of a
// particular type. The template argument \a MT represents the type of the matrix operand.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
class OperationTest
{
 private:
   //**Type definitions****************************************************************************
   using ET = blaze::ElementType_t<MT>;  //!< Element type.

   using OMT  = blaze::OppositeType_t<MT>;    //!< Matrix type with opposite storage order
   using TMT  = blaze::TransposeType_t<MT>;   //!< Transpose matrix type.
   using TOMT = blaze::TransposeType_t<OMT>;  //!< Transpose matrix type with opposite storage order

   //! Sparse result type
   using SRE = blaze::RepeatTrait_t<MT,R0,R1>;

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

   //! Reference type.
   using RT = blaze::DynamicMatrix<ET,false>;

   //! Reference result type
   using RRE = MatchSymmetry_t< SRE, blaze::RepeatTrait_t<RT,R0,R1> >;

   //! Type of the matrix repeater expression (runtime argument)
   using MatRepeatExprType1 =
      blaze::RemoveCVRef_t< decltype( blaze::repeat( std::declval<MT>(), R0, R1 ) ) >;

   //! Type of the matrix repeater expression (compile time argument)
   using MatRepeatExprType2 =
      blaze::RemoveCVRef_t< decltype( blaze::repeat<R0,R1>( std::declval<MT>() ) ) >;

   //! Type of the transpose matrix repeater expression (runtime argument)
   using TMatRepeatExprType1 =
      blaze::RemoveCVRef_t< decltype( blaze::repeat( std::declval<OMT>(), R0, R1 ) ) >;

   //! Type of the transpose matrix repeater expression (compile time argument)
   using TMatRepeatExprType2 =
      blaze::RemoveCVRef_t< decltype( blaze::repeat<R0,R1>( std::declval<OMT>() ) ) >;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit OperationTest( const Creator<MT>& creator );
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
   template< typename T > void checkResults();
   template< typename T > void checkTransposeResults();
   void checkExceptionMessage( const std::exception& ex, const std::string& message );
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
   MT    mat_;      //!< The sparse matrix operand.
   OMT   omat_;     //!< The sparse matrix with opposite storage order.
   DRE   dres_;     //!< The dense result matrix.
   SRE   sres_;     //!< The sparse result matrix.
   ODRE  odres_;    //!< The dense result matrix with opposite storage order.
   OSRE  osres_;    //!< The sparse result matrix with opposite storage order.
   TDRE  tdres_;    //!< The transpose dense result matrix.
   TSRE  tsres_;    //!< The transpose sparse result matrix.
   TODRE todres_;   //!< The transpose dense result matrix with opposite storage order.
   TOSRE tosres_;   //!< The transpose sparse result matrix with opposite storage order.
   RT    refmat_;   //!< The reference matrix.
   RRE   refres_;   //!< The reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT    );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OMT   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TMT   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TOMT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( RT    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( RRE   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( ODRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OSRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TODRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TOSRE );

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( MT    );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OMT   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TMT   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TOMT  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( RT    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( RRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ODRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OSRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TODRE );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TOSRE );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DRE   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ODRE  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TDRE  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TODRE );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SRE   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OSRE  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TSRE  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TOSRE );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RRE   );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET , blaze::ElementType_t<OMT>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET , blaze::ElementType_t<TMT>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET , blaze::ElementType_t<TOMT>  );
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
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT , blaze::OppositeType_t<OMT>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT , blaze::TransposeType_t<TMT>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DRE, blaze::OppositeType_t<ODRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DRE, blaze::TransposeType_t<TDRE> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SRE, blaze::OppositeType_t<OSRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SRE, blaze::TransposeType_t<TSRE> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( MatRepeatExprType1, blaze::ResultType_t<MatRepeatExprType1>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MatRepeatExprType1, blaze::OppositeType_t<MatRepeatExprType1>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MatRepeatExprType1, blaze::TransposeType_t<MatRepeatExprType1> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( MatRepeatExprType2, blaze::ResultType_t<MatRepeatExprType2>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MatRepeatExprType2, blaze::OppositeType_t<MatRepeatExprType2>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( MatRepeatExprType2, blaze::TransposeType_t<MatRepeatExprType2> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( TMatRepeatExprType1, blaze::ResultType_t<TMatRepeatExprType1>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TMatRepeatExprType1, blaze::OppositeType_t<TMatRepeatExprType1>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TMatRepeatExprType1, blaze::TransposeType_t<TMatRepeatExprType1> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( TMatRepeatExprType2, blaze::ResultType_t<TMatRepeatExprType2>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TMatRepeatExprType2, blaze::OppositeType_t<TMatRepeatExprType2>  );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TMatRepeatExprType2, blaze::TransposeType_t<TMatRepeatExprType2> );

   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( MatRepeatExprType1 , blaze::BaseType_t<MatRepeatExprType1 > );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( MatRepeatExprType2 , blaze::BaseType_t<MatRepeatExprType2 > );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TMatRepeatExprType1, blaze::BaseType_t<TMatRepeatExprType1> );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TMatRepeatExprType2, blaze::BaseType_t<TMatRepeatExprType2> );
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
/*!\brief Constructor for the sparse matrix repeat operation test.
//
// \param creator The creator for sparse matrix operand.
// \exception std::runtime_error Operation error detected.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
OperationTest<MT,R0,R1>::OperationTest( const Creator<MT>& creator )
   : mat_( creator() )  // The sparse matrix operand
   , omat_( mat_ )      // The sparse matrix with opposite storage order
   , dres_()            // The dense result matrix
   , sres_()            // The sparse result matrix
   , odres_()           // The dense result matrix with opposite storage order
   , osres_()           // The sparse result matrix with opposite storage order
   , tdres_()           // The transpose dense result matrix
   , tsres_()           // The transpose sparse result matrix
   , todres_()          // The transpose dense result matrix with opposite storage order
   , tosres_()          // The transpose sparse result matrix with opposite storage order
   , refmat_( mat_ )    // The reference matrix
   , refres_()          // The reference result
   , test_()            // Label of the currently performed test
   , error_()           // Description of the current error type
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
/*!\brief Tests on the initial status of the matrix.
//
// \return void
// \exception std::runtime_error Initialization error detected.
//
// This function runs tests on the initial status of the matrix. In case any initialization
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testInitialStatus()
{
   //=====================================================================================
   // Performing initial tests with the row-major types
   //=====================================================================================

   // Checking the number of rows of the matrix operand
   if( mat_.rows() != refmat_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of sparse matrix operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Detected number of rows = " << mat_.rows() << "\n"
          << "   Expected number of rows = " << refmat_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the matrix operand
   if( mat_.columns() != refmat_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of sparse matrix operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Detected number of columns = " << mat_.columns() << "\n"
          << "   Expected number of columns = " << refmat_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the matrix operand
   if( !isEqual( mat_, refmat_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of sparse matrix operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Current initialization:\n" << mat_ << "\n"
          << "   Expected initialization:\n" << refmat_ << "\n";
      throw std::runtime_error( oss.str() );
   }


   //=====================================================================================
   // Performing initial tests with the column-major types
   //=====================================================================================

   // Checking the number of rows of the matrix operand
   if( omat_.rows() != refmat_.rows() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of sparse matrix operand\n"
          << " Error: Invalid number of rows\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse matrix type:\n"
          << "     " << typeid( OMT ).name() << "\n"
          << "   Detected number of rows = " << omat_.rows() << "\n"
          << "   Expected number of rows = " << refmat_.rows() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the number of columns of the matrix operand
   if( omat_.columns() != refmat_.columns() ) {
      std::ostringstream oss;
      oss << " Test: Initial size comparison of sparse matrix operand\n"
          << " Error: Invalid number of columns\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse matrix type:\n"
          << "     " << typeid( OMT ).name() << "\n"
          << "   Detected number of columns = " << omat_.columns() << "\n"
          << "   Expected number of columns = " << refmat_.columns() << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Checking the initialization of the matrix operand
   if( !isEqual( omat_, refmat_ ) ) {
      std::ostringstream oss;
      oss << " Test: Initial test of initialization of sparse matrix operand\n"
          << " Error: Invalid matrix initialization\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse matrix type:\n"
          << "     " << typeid( OMT ).name() << "\n"
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
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testAssignment()
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
      oss << " Test: Checking the assignment result of row-major sparse matrix operand\n"
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
   // Performing an assignment with the transpose type
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

   if( !isEqual( mat_, refmat_ ) ) {
      std::ostringstream oss;
      oss << " Test: Checking the assignment result of column-major sparse matrix operand\n"
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
/*!\brief Testing the explicit evaluation.
//
// \return void
// \exception std::runtime_error Evaluation error detected.
//
// This function tests the explicit evaluation. In case any error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testEvaluation()
{
   using blaze::IsRowMajorMatrix;
   using blaze::repeat;


   //=====================================================================================
   // Testing the evaluation with a row-major matrix
   //=====================================================================================

   {
      const auto res   ( evaluate( repeat( mat_, R0, R1 ) ) );
      const auto refres( evaluate( repeat( refmat_, R0, R1 ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrix (runtime)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse row-major matrix type:\n"
             << "     " << typeid( mat_ ).name() << "\n"
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
      const auto res   ( evaluate( repeat<R0,R1>( mat_ ) ) );
      const auto refres( evaluate( repeat<R0,R1>( refmat_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrix (compile time)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse row-major matrix type:\n"
             << "     " << typeid( mat_ ).name() << "\n"
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
      const auto res   ( evaluate( repeat( eval( mat_ ), R0, R1 ) ) );
      const auto refres( evaluate( repeat( eval( refmat_ ), R0, R1 ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated matrix (runtime)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse row-major matrix type:\n"
             << "     " << typeid( mat_ ).name() << "\n"
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
      const auto res   ( evaluate( repeat<R0,R1>( eval( mat_ ) ) ) );
      const auto refres( evaluate( repeat<R0,R1>( eval( refmat_ ) ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated matrix (compile time)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse row-major matrix type:\n"
             << "     " << typeid( mat_ ).name() << "\n"
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
   // Testing the evaluation with a column-major matrix
   //=====================================================================================

   {
      const auto res   ( evaluate( repeat( omat_, R0, R1 ) ) );
      const auto refres( evaluate( repeat( refmat_, R0, R1 ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrix (runtime)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse column-major matrix type:\n"
             << "     " << typeid( omat_ ).name() << "\n"
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
      const auto res   ( evaluate( repeat<R0,R1>( omat_ ) ) );
      const auto refres( evaluate( repeat<R0,R1>( refmat_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given matrix (compile time)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse column-major matrix type:\n"
             << "     " << typeid( omat_ ).name() << "\n"
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
      const auto res   ( evaluate( repeat( eval( omat_ ), R0, R1 ) ) );
      const auto refres( evaluate( repeat( eval( refmat_ ), R0, R1 ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated matrix (runtime)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse column-major matrix type:\n"
             << "     " << typeid( omat_ ).name() << "\n"
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
      const auto res   ( evaluate( repeat<R0,R1>( eval( omat_ ) ) ) );
      const auto refres( evaluate( repeat<R0,R1>( eval( refmat_ ) ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated matrix (compile time)\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse column-major matrix type:\n"
             << "     " << typeid( omat_ ).name() << "\n"
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
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testElementAccess()
{
   using blaze::equal;
   using blaze::repeat;


   //=====================================================================================
   // Testing the element access with a row-major matrix
   //=====================================================================================

   if( mat_.rows() > 0UL && mat_.columns() > 0UL )
   {
      const size_t m( mat_.rows()*R0    - 1UL );
      const size_t n( mat_.columns()*R1 - 1UL );

      if( !equal( repeat( mat_, R0, R1 )(m,n), repeat( refmat_, R0, R1 )(m,n) ) ||
          !equal( repeat( mat_, R0, R1 ).at(m,n), repeat( refmat_, R0, R1 ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of repeater expression (runtime)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse row-major matrix type:\n"
             << "     " << typeid( MT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( repeat<R0,R1>( mat_ )(m,n), repeat<R0,R1>( refmat_ )(m,n) ) ||
          !equal( repeat<R0,R1>( mat_ ).at(m,n), repeat<R0,R1>( refmat_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of repeater expression (compile time)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse row-major matrix type:\n"
             << "     " << typeid( MT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( repeat( eval( mat_ ), R0, R1 )(m,n), repeat( eval( refmat_ ), R0, R1 )(m,n) ) ||
          !equal( repeat( eval( mat_ ), R0, R1 ).at(m,n), repeat( eval( refmat_ ), R0, R1 ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of evaluated repeater expression (runtime)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse row-major matrix type:\n"
             << "     " << typeid( MT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( repeat<R0,R1>( eval( mat_ ) )(m,n), repeat<R0,R1>( eval( refmat_ ) )(m,n) ) ||
          !equal( repeat<R0,R1>( eval( mat_ ) ).at(m,n), repeat<R0,R1>( eval( refmat_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of evaluated repeater expression (compile time)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse row-major matrix type:\n"
             << "     " << typeid( MT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      repeat( mat_, R0, R1 ).at( 0UL, mat_.columns()*R1 );

      std::ostringstream oss;
      oss << " Test : Checked element access of repeater expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse row-major matrix type:\n"
          << "     " << typeid( MT ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      repeat<R0,R1>( mat_ ).at( mat_.rows()*R0, 0UL );

      std::ostringstream oss;
      oss << " Test : Checked element access of repeater expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse row-major matrix type:\n"
          << "     " << typeid( MT ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}


   //=====================================================================================
   // Testing the element access with a column-major matrix
   //=====================================================================================

   if( omat_.rows() > 0UL && omat_.columns() > 0UL )
   {
      const size_t m( omat_.rows()*R0    - 1UL );
      const size_t n( omat_.columns()*R1 - 1UL );

      if( !equal( repeat( omat_, R0, R1 )(m,n), repeat( refmat_, R0, R1 )(m,n) ) ||
          !equal( repeat( omat_, R0, R1 ).at(m,n), repeat( refmat_, R0, R1 ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of repeater expression (runtime)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse column-major matrix type:\n"
             << "     " << typeid( OMT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( repeat<R0,R1>( omat_ )(m,n), repeat<R0,R1>( refmat_ )(m,n) ) ||
          !equal( repeat<R0,R1>( omat_ ).at(m,n), repeat<R0,R1>( refmat_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of repeater expression (compile time)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse column-major matrix type:\n"
             << "     " << typeid( OMT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( repeat( eval( omat_ ), R0, R1 )(m,n), repeat( eval( refmat_ ), R0, R1 )(m,n) ) ||
          !equal( repeat( eval( omat_ ), R0, R1 ).at(m,n), repeat( eval( refmat_ ), R0, R1 ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of evaluated repeater expression (runtime)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse column-major matrix type:\n"
             << "     " << typeid( OMT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( repeat<R0,R1>( eval( omat_ ) )(m,n), repeat<R0,R1>( eval( refmat_ ) )(m,n) ) ||
          !equal( repeat<R0,R1>( eval( omat_ ) ).at(m,n), repeat<R0,R1>( eval( refmat_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of evaluated repeater expression (compile time)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Sparse column-major matrix type:\n"
             << "     " << typeid( OMT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      repeat( omat_, R0, R1 ).at( 0UL, omat_.columns()*R1 );

      std::ostringstream oss;
      oss << " Test : Checked element access of repeater expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse column-major matrix type:\n"
          << "     " << typeid( OMT ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   try {
      repeat<R0,R1>( omat_ ).at( omat_.rows()*R0, 0UL );

      std::ostringstream oss;
      oss << " Test : Checked element access of repeater expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse column-major matrix type:\n"
          << "     " << typeid( OMT ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the matrix repeat operation with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case
// any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      using blaze::repeat;


      //=====================================================================================
      // Repeat operation
      //=====================================================================================

      // Repeat operation with the given matrix (runtime)
      {
         test_  = "Repeat operation with the given matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat( mat_, R0, R1 );
            odres_  = repeat( mat_, R0, R1 );
            sres_   = repeat( mat_, R0, R1 );
            osres_  = repeat( mat_, R0, R1 );
            refres_ = repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = repeat( omat_, R0, R1 );
            odres_  = repeat( omat_, R0, R1 );
            sres_   = repeat( omat_, R0, R1 );
            osres_  = repeat( omat_, R0, R1 );
            refres_ = repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Repeat operation with the given matrix (compile time)
      {
         test_  = "Repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat<R0,R1>( mat_ );
            odres_  = repeat<R0,R1>( mat_ );
            sres_   = repeat<R0,R1>( mat_ );
            osres_  = repeat<R0,R1>( mat_ );
            refres_ = repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = repeat<R0,R1>( omat_ );
            odres_  = repeat<R0,R1>( omat_ );
            sres_   = repeat<R0,R1>( omat_ );
            osres_  = repeat<R0,R1>( omat_ );
            refres_ = repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Repeat operation with evaluated matrix (runtime)
      {
         test_  = "Repeat operation with evaluated matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat( eval( mat_ ), R0, R1 );
            odres_  = repeat( eval( mat_ ), R0, R1 );
            sres_   = repeat( eval( mat_ ), R0, R1 );
            osres_  = repeat( eval( mat_ ), R0, R1 );
            refres_ = repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = repeat( eval( omat_ ), R0, R1 );
            odres_  = repeat( eval( omat_ ), R0, R1 );
            sres_   = repeat( eval( omat_ ), R0, R1 );
            osres_  = repeat( eval( omat_ ), R0, R1 );
            refres_ = repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Repeat operation with evaluated matrix (compile time)
      {
         test_  = "Repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat<R0,R1>( eval( mat_ ) );
            odres_  = repeat<R0,R1>( eval( mat_ ) );
            sres_   = repeat<R0,R1>( eval( mat_ ) );
            osres_  = repeat<R0,R1>( eval( mat_ ) );
            refres_ = repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = repeat<R0,R1>( eval( omat_ ) );
            odres_  = repeat<R0,R1>( eval( omat_ ) );
            sres_   = repeat<R0,R1>( eval( omat_ ) );
            osres_  = repeat<R0,R1>( eval( omat_ ) );
            refres_ = repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Repeat with addition assignment
      //=====================================================================================

      // Repeat with addition assignment with the given matrix (runtime)
      {
         test_  = "Repeat with addition assignment with the given matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat( mat_, R0, R1 );
            odres_  += repeat( mat_, R0, R1 );
            sres_   += repeat( mat_, R0, R1 );
            osres_  += repeat( mat_, R0, R1 );
            refres_ += repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += repeat( omat_, R0, R1 );
            odres_  += repeat( omat_, R0, R1 );
            sres_   += repeat( omat_, R0, R1 );
            osres_  += repeat( omat_, R0, R1 );
            refres_ += repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Repeat with addition assignment with the given matrix (compile time)
      {
         test_  = "Repeat with addition assignment with the given matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat<R0,R1>( mat_ );
            odres_  += repeat<R0,R1>( mat_ );
            sres_   += repeat<R0,R1>( mat_ );
            osres_  += repeat<R0,R1>( mat_ );
            refres_ += repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += repeat<R0,R1>( omat_ );
            odres_  += repeat<R0,R1>( omat_ );
            sres_   += repeat<R0,R1>( omat_ );
            osres_  += repeat<R0,R1>( omat_ );
            refres_ += repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Repeat with addition assignment with evaluated matrix (runtime)
      {
         test_  = "Repeat with addition assignment with evaluated matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat( eval( mat_ ), R0, R1 );
            odres_  += repeat( eval( mat_ ), R0, R1 );
            sres_   += repeat( eval( mat_ ), R0, R1 );
            osres_  += repeat( eval( mat_ ), R0, R1 );
            refres_ += repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += repeat( eval( omat_ ), R0, R1 );
            odres_  += repeat( eval( omat_ ), R0, R1 );
            sres_   += repeat( eval( omat_ ), R0, R1 );
            osres_  += repeat( eval( omat_ ), R0, R1 );
            refres_ += repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Repeat with addition assignment with evaluated matrix (compile time)
      {
         test_  = "Repeat with addition assignment with the given matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat<R0,R1>( eval( mat_ ) );
            odres_  += repeat<R0,R1>( eval( mat_ ) );
            sres_   += repeat<R0,R1>( eval( mat_ ) );
            osres_  += repeat<R0,R1>( eval( mat_ ) );
            refres_ += repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += repeat<R0,R1>( eval( omat_ ) );
            odres_  += repeat<R0,R1>( eval( omat_ ) );
            sres_   += repeat<R0,R1>( eval( omat_ ) );
            osres_  += repeat<R0,R1>( eval( omat_ ) );
            refres_ += repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Repeat with subtraction assignment
      //=====================================================================================

      // Repeat with subtraction assignment with the given matrix (runtime)
      {
         test_  = "Repeat with subtraction assignment with the given matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat( mat_, R0, R1 );
            odres_  -= repeat( mat_, R0, R1 );
            sres_   -= repeat( mat_, R0, R1 );
            osres_  -= repeat( mat_, R0, R1 );
            refres_ -= repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= repeat( omat_, R0, R1 );
            odres_  -= repeat( omat_, R0, R1 );
            sres_   -= repeat( omat_, R0, R1 );
            osres_  -= repeat( omat_, R0, R1 );
            refres_ -= repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Repeat with subtraction assignment with the given matrix (compile time)
      {
         test_  = "Repeat with subtraction assignment with the given matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat<R0,R1>( mat_ );
            odres_  -= repeat<R0,R1>( mat_ );
            sres_   -= repeat<R0,R1>( mat_ );
            osres_  -= repeat<R0,R1>( mat_ );
            refres_ -= repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= repeat<R0,R1>( omat_ );
            odres_  -= repeat<R0,R1>( omat_ );
            sres_   -= repeat<R0,R1>( omat_ );
            osres_  -= repeat<R0,R1>( omat_ );
            refres_ -= repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Repeat with subtraction assignment with evaluated matrix (runtime)
      {
         test_  = "Repeat with subtraction assignment with evaluated matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat( eval( mat_ ), R0, R1 );
            odres_  -= repeat( eval( mat_ ), R0, R1 );
            sres_   -= repeat( eval( mat_ ), R0, R1 );
            osres_  -= repeat( eval( mat_ ), R0, R1 );
            refres_ -= repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= repeat( eval( omat_ ), R0, R1 );
            odres_  -= repeat( eval( omat_ ), R0, R1 );
            sres_   -= repeat( eval( omat_ ), R0, R1 );
            osres_  -= repeat( eval( omat_ ), R0, R1 );
            refres_ -= repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Repeat with subtraction assignment with evaluated matrix (compile time)
      {
         test_  = "Repeat with subtraction assignment with the given matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat<R0,R1>( eval( mat_ ) );
            odres_  -= repeat<R0,R1>( eval( mat_ ) );
            sres_   -= repeat<R0,R1>( eval( mat_ ) );
            osres_  -= repeat<R0,R1>( eval( mat_ ) );
            refres_ -= repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= repeat<R0,R1>( eval( omat_ ) );
            odres_  -= repeat<R0,R1>( eval( omat_ ) );
            sres_   -= repeat<R0,R1>( eval( omat_ ) );
            osres_  -= repeat<R0,R1>( eval( omat_ ) );
            refres_ -= repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Repeat with Schur product assignment
      //=====================================================================================

      // Repeat with Schur product assignment with the given matrix (runtime)
      {
         test_  = "Repeat with Schur product assignment with the given matrix (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= repeat( mat_, R0, R1 );
            odres_  %= repeat( mat_, R0, R1 );
            sres_   %= repeat( mat_, R0, R1 );
            osres_  %= repeat( mat_, R0, R1 );
            refres_ %= repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= repeat( omat_, R0, R1 );
            odres_  %= repeat( omat_, R0, R1 );
            sres_   %= repeat( omat_, R0, R1 );
            osres_  %= repeat( omat_, R0, R1 );
            refres_ %= repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Repeat with Schur product assignment with the given matrix (compile time)
      {
         test_  = "Repeat with Schur product assignment with the given matrix (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= repeat<R0,R1>( mat_ );
            odres_  %= repeat<R0,R1>( mat_ );
            sres_   %= repeat<R0,R1>( mat_ );
            osres_  %= repeat<R0,R1>( mat_ );
            refres_ %= repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= repeat<R0,R1>( omat_ );
            odres_  %= repeat<R0,R1>( omat_ );
            sres_   %= repeat<R0,R1>( omat_ );
            osres_  %= repeat<R0,R1>( omat_ );
            refres_ %= repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Repeat with Schur product assignment with evaluated matrix (runtime)
      {
         test_  = "Repeat with Schur product assignment with evaluated matrix (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= repeat( eval( mat_ ), R0, R1 );
            odres_  %= repeat( eval( mat_ ), R0, R1 );
            sres_   %= repeat( eval( mat_ ), R0, R1 );
            osres_  %= repeat( eval( mat_ ), R0, R1 );
            refres_ %= repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= repeat( eval( omat_ ), R0, R1 );
            odres_  %= repeat( eval( omat_ ), R0, R1 );
            sres_   %= repeat( eval( omat_ ), R0, R1 );
            osres_  %= repeat( eval( omat_ ), R0, R1 );
            refres_ %= repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Repeat with Schur product assignment with evaluated matrix (compile time)
      {
         test_  = "Repeat with Schur product assignment with the given matrix (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= repeat<R0,R1>( eval( mat_ ) );
            odres_  %= repeat<R0,R1>( eval( mat_ ) );
            sres_   %= repeat<R0,R1>( eval( mat_ ) );
            osres_  %= repeat<R0,R1>( eval( mat_ ) );
            refres_ %= repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= repeat<R0,R1>( eval( omat_ ) );
            odres_  %= repeat<R0,R1>( eval( omat_ ) );
            sres_   %= repeat<R0,R1>( eval( omat_ ) );
            osres_  %= repeat<R0,R1>( eval( omat_ ) );
            refres_ %= repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the negated sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the negated matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment, and division assignment.
// In case any error resulting from the repeat operation or the subsequent assignment is
// detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      using blaze::repeat;


      //=====================================================================================
      // Negated repeat operation
      //=====================================================================================

      // Negated repeat operation with the given matrix (runtime)
      {
         test_  = "Negated repeat operation with the given matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = -repeat( mat_, R0, R1 );
            odres_  = -repeat( mat_, R0, R1 );
            sres_   = -repeat( mat_, R0, R1 );
            osres_  = -repeat( mat_, R0, R1 );
            refres_ = -repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = -repeat( omat_, R0, R1 );
            odres_  = -repeat( omat_, R0, R1 );
            sres_   = -repeat( omat_, R0, R1 );
            osres_  = -repeat( omat_, R0, R1 );
            refres_ = -repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Negated repeat operation with the given matrix (compile time)
      {
         test_  = "Negated repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = -repeat<R0,R1>( mat_ );
            odres_  = -repeat<R0,R1>( mat_ );
            sres_   = -repeat<R0,R1>( mat_ );
            osres_  = -repeat<R0,R1>( mat_ );
            refres_ = -repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = -repeat<R0,R1>( omat_ );
            odres_  = -repeat<R0,R1>( omat_ );
            sres_   = -repeat<R0,R1>( omat_ );
            osres_  = -repeat<R0,R1>( omat_ );
            refres_ = -repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Negated repeat operation with evaluated matrix (runtime)
      {
         test_  = "Negated repeat operation with evaluated matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = -repeat( eval( mat_ ), R0, R1 );
            odres_  = -repeat( eval( mat_ ), R0, R1 );
            sres_   = -repeat( eval( mat_ ), R0, R1 );
            osres_  = -repeat( eval( mat_ ), R0, R1 );
            refres_ = -repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = -repeat( eval( omat_ ), R0, R1 );
            odres_  = -repeat( eval( omat_ ), R0, R1 );
            sres_   = -repeat( eval( omat_ ), R0, R1 );
            osres_  = -repeat( eval( omat_ ), R0, R1 );
            refres_ = -repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Negated repeat operation with evaluated matrix (compile time)
      {
         test_  = "Negated repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = -repeat<R0,R1>( eval( mat_ ) );
            odres_  = -repeat<R0,R1>( eval( mat_ ) );
            sres_   = -repeat<R0,R1>( eval( mat_ ) );
            osres_  = -repeat<R0,R1>( eval( mat_ ) );
            refres_ = -repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = -repeat<R0,R1>( eval( omat_ ) );
            odres_  = -repeat<R0,R1>( eval( omat_ ) );
            sres_   = -repeat<R0,R1>( eval( omat_ ) );
            osres_  = -repeat<R0,R1>( eval( omat_ ) );
            refres_ = -repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Negated repeat with addition assignment
      //=====================================================================================

      // Negated repeat with addition assignment with the given matrix (runtime)
      {
         test_  = "Negated repeat with addition assignment with the given matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += -repeat( mat_, R0, R1 );
            odres_  += -repeat( mat_, R0, R1 );
            sres_   += -repeat( mat_, R0, R1 );
            osres_  += -repeat( mat_, R0, R1 );
            refres_ += -repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += -repeat( omat_, R0, R1 );
            odres_  += -repeat( omat_, R0, R1 );
            sres_   += -repeat( omat_, R0, R1 );
            osres_  += -repeat( omat_, R0, R1 );
            refres_ += -repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Negated repeat with addition assignment with the given matrix (compile time)
      {
         test_  = "Negated repeat with addition assignment with the given matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += -repeat<R0,R1>( mat_ );
            odres_  += -repeat<R0,R1>( mat_ );
            sres_   += -repeat<R0,R1>( mat_ );
            osres_  += -repeat<R0,R1>( mat_ );
            refres_ += -repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += -repeat<R0,R1>( omat_ );
            odres_  += -repeat<R0,R1>( omat_ );
            sres_   += -repeat<R0,R1>( omat_ );
            osres_  += -repeat<R0,R1>( omat_ );
            refres_ += -repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Negated repeat with addition assignment with evaluated matrix (runtime)
      {
         test_  = "Negated repeat with addition assignment with evaluated matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += -repeat( eval( mat_ ), R0, R1 );
            odres_  += -repeat( eval( mat_ ), R0, R1 );
            sres_   += -repeat( eval( mat_ ), R0, R1 );
            osres_  += -repeat( eval( mat_ ), R0, R1 );
            refres_ += -repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += -repeat( eval( omat_ ), R0, R1 );
            odres_  += -repeat( eval( omat_ ), R0, R1 );
            sres_   += -repeat( eval( omat_ ), R0, R1 );
            osres_  += -repeat( eval( omat_ ), R0, R1 );
            refres_ += -repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Negated repeat with addition assignment with evaluated matrix (compile time)
      {
         test_  = "Negated repeat with addition assignment with the given matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += -repeat<R0,R1>( eval( mat_ ) );
            odres_  += -repeat<R0,R1>( eval( mat_ ) );
            sres_   += -repeat<R0,R1>( eval( mat_ ) );
            osres_  += -repeat<R0,R1>( eval( mat_ ) );
            refres_ += -repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += -repeat<R0,R1>( eval( omat_ ) );
            odres_  += -repeat<R0,R1>( eval( omat_ ) );
            sres_   += -repeat<R0,R1>( eval( omat_ ) );
            osres_  += -repeat<R0,R1>( eval( omat_ ) );
            refres_ += -repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Negated repeat with subtraction assignment
      //=====================================================================================

      // Negated repeat with subtraction assignment with the given matrix (runtime)
      {
         test_  = "Negated repeat with subtraction assignment with the given matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= -repeat( mat_, R0, R1 );
            odres_  -= -repeat( mat_, R0, R1 );
            sres_   -= -repeat( mat_, R0, R1 );
            osres_  -= -repeat( mat_, R0, R1 );
            refres_ -= -repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= -repeat( omat_, R0, R1 );
            odres_  -= -repeat( omat_, R0, R1 );
            sres_   -= -repeat( omat_, R0, R1 );
            osres_  -= -repeat( omat_, R0, R1 );
            refres_ -= -repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Negated repeat with subtraction assignment with the given matrix (compile time)
      {
         test_  = "Negated repeat with subtraction assignment with the given matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= -repeat<R0,R1>( mat_ );
            odres_  -= -repeat<R0,R1>( mat_ );
            sres_   -= -repeat<R0,R1>( mat_ );
            osres_  -= -repeat<R0,R1>( mat_ );
            refres_ -= -repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= -repeat<R0,R1>( omat_ );
            odres_  -= -repeat<R0,R1>( omat_ );
            sres_   -= -repeat<R0,R1>( omat_ );
            osres_  -= -repeat<R0,R1>( omat_ );
            refres_ -= -repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Negated repeat with subtraction assignment with evaluated matrix (runtime)
      {
         test_  = "Negated repeat with subtraction assignment with evaluated matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= -repeat( eval( mat_ ), R0, R1 );
            odres_  -= -repeat( eval( mat_ ), R0, R1 );
            sres_   -= -repeat( eval( mat_ ), R0, R1 );
            osres_  -= -repeat( eval( mat_ ), R0, R1 );
            refres_ -= -repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= -repeat( eval( omat_ ), R0, R1 );
            odres_  -= -repeat( eval( omat_ ), R0, R1 );
            sres_   -= -repeat( eval( omat_ ), R0, R1 );
            osres_  -= -repeat( eval( omat_ ), R0, R1 );
            refres_ -= -repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Negated repeat with subtraction assignment with evaluated matrix (compile time)
      {
         test_  = "Negated repeat with subtraction assignment with the given matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= -repeat<R0,R1>( eval( mat_ ) );
            odres_  -= -repeat<R0,R1>( eval( mat_ ) );
            sres_   -= -repeat<R0,R1>( eval( mat_ ) );
            osres_  -= -repeat<R0,R1>( eval( mat_ ) );
            refres_ -= -repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= -repeat<R0,R1>( eval( omat_ ) );
            odres_  -= -repeat<R0,R1>( eval( omat_ ) );
            sres_   -= -repeat<R0,R1>( eval( omat_ ) );
            osres_  -= -repeat<R0,R1>( eval( omat_ ) );
            refres_ -= -repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Negated repeat with Schur product assignment
      //=====================================================================================

      // Negated repeat with Schur product assignment with the given matrix (runtime)
      {
         test_  = "Negated repeat with Schur product assignment with the given matrix (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= -repeat( mat_, R0, R1 );
            odres_  %= -repeat( mat_, R0, R1 );
            sres_   %= -repeat( mat_, R0, R1 );
            osres_  %= -repeat( mat_, R0, R1 );
            refres_ %= -repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= -repeat( omat_, R0, R1 );
            odres_  %= -repeat( omat_, R0, R1 );
            sres_   %= -repeat( omat_, R0, R1 );
            osres_  %= -repeat( omat_, R0, R1 );
            refres_ %= -repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Negated repeat with Schur product assignment with the given matrix (compile time)
      {
         test_  = "Negated repeat with Schur product assignment with the given matrix (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= -repeat<R0,R1>( mat_ );
            odres_  %= -repeat<R0,R1>( mat_ );
            sres_   %= -repeat<R0,R1>( mat_ );
            osres_  %= -repeat<R0,R1>( mat_ );
            refres_ %= -repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= -repeat<R0,R1>( omat_ );
            odres_  %= -repeat<R0,R1>( omat_ );
            sres_   %= -repeat<R0,R1>( omat_ );
            osres_  %= -repeat<R0,R1>( omat_ );
            refres_ %= -repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Negated repeat with Schur product assignment with evaluated matrix (runtime)
      {
         test_  = "Negated repeat with Schur product assignment with evaluated matrix (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= -repeat( eval( mat_ ), R0, R1 );
            odres_  %= -repeat( eval( mat_ ), R0, R1 );
            sres_   %= -repeat( eval( mat_ ), R0, R1 );
            osres_  %= -repeat( eval( mat_ ), R0, R1 );
            refres_ %= -repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= -repeat( eval( omat_ ), R0, R1 );
            odres_  %= -repeat( eval( omat_ ), R0, R1 );
            sres_   %= -repeat( eval( omat_ ), R0, R1 );
            osres_  %= -repeat( eval( omat_ ), R0, R1 );
            refres_ %= -repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Negated repeat with Schur product assignment with evaluated matrix (compile time)
      {
         test_  = "Negated repeat with Schur product assignment with the given matrix (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= -repeat<R0,R1>( eval( mat_ ) );
            odres_  %= -repeat<R0,R1>( eval( mat_ ) );
            sres_   %= -repeat<R0,R1>( eval( mat_ ) );
            osres_  %= -repeat<R0,R1>( eval( mat_ ) );
            refres_ %= -repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= -repeat<R0,R1>( eval( omat_ ) );
            odres_  %= -repeat<R0,R1>( eval( omat_ ) );
            sres_   %= -repeat<R0,R1>( eval( omat_ ) );
            osres_  %= -repeat<R0,R1>( eval( omat_ ) );
            refres_ %= -repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the scaled sparse matrix repeat operation.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the scaled matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT   // Type of the sparse matrix
        , size_t R0     // Compile time row-wise repetitions
        , size_t R1 >   // Compile time column-wise repetitions
template< typename T >  // Type of the scalar
void OperationTest<MT,R0,R1>::testScaledOperation( T scalar )
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

      // Scaled repeat operation with the given matrix (s*OP, runtime)
      {
         test_  = "Scaled repeat operation with the given matrix (s*OP, runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = scalar * repeat( mat_, R0, R1 );
            odres_  = scalar * repeat( mat_, R0, R1 );
            sres_   = scalar * repeat( mat_, R0, R1 );
            osres_  = scalar * repeat( mat_, R0, R1 );
            refres_ = scalar * repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = scalar * repeat( omat_, R0, R1 );
            odres_  = scalar * repeat( omat_, R0, R1 );
            sres_   = scalar * repeat( omat_, R0, R1 );
            osres_  = scalar * repeat( omat_, R0, R1 );
            refres_ = scalar * repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat operation with the given matrix (s*OP, compile time)
      {
         test_  = "Scaled repeat operation with the given matrix (s*OP, compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = scalar * repeat<R0,R1>( mat_ );
            odres_  = scalar * repeat<R0,R1>( mat_ );
            sres_   = scalar * repeat<R0,R1>( mat_ );
            osres_  = scalar * repeat<R0,R1>( mat_ );
            refres_ = scalar * repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = scalar * repeat<R0,R1>( omat_ );
            odres_  = scalar * repeat<R0,R1>( omat_ );
            sres_   = scalar * repeat<R0,R1>( omat_ );
            osres_  = scalar * repeat<R0,R1>( omat_ );
            refres_ = scalar * repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat operation with evaluated matrix (s*OP, runtime)
      {
         test_  = "Scaled repeat operation with evaluated matrix (s*OP, runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = scalar * repeat( eval( mat_ ), R0, R1 );
            odres_  = scalar * repeat( eval( mat_ ), R0, R1 );
            sres_   = scalar * repeat( eval( mat_ ), R0, R1 );
            osres_  = scalar * repeat( eval( mat_ ), R0, R1 );
            refres_ = scalar * repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = scalar * repeat( eval( omat_ ), R0, R1 );
            odres_  = scalar * repeat( eval( omat_ ), R0, R1 );
            sres_   = scalar * repeat( eval( omat_ ), R0, R1 );
            osres_  = scalar * repeat( eval( omat_ ), R0, R1 );
            refres_ = scalar * repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat operation with evaluated matrix (s*OP, compile time)
      {
         test_  = "Scaled repeat operation with the given matrix (s*OP, compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = scalar * repeat<R0,R1>( eval( mat_ ) );
            odres_  = scalar * repeat<R0,R1>( eval( mat_ ) );
            sres_   = scalar * repeat<R0,R1>( eval( mat_ ) );
            osres_  = scalar * repeat<R0,R1>( eval( mat_ ) );
            refres_ = scalar * repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = scalar * repeat<R0,R1>( eval( omat_ ) );
            odres_  = scalar * repeat<R0,R1>( eval( omat_ ) );
            sres_   = scalar * repeat<R0,R1>( eval( omat_ ) );
            osres_  = scalar * repeat<R0,R1>( eval( omat_ ) );
            refres_ = scalar * repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Scaled repeat operation (OP*s)
      //=====================================================================================

      // Scaled repeat operation with the given matrix (OP*s, runtime)
      {
         test_  = "Scaled repeat operation with the given matrix (OP*s, runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat( mat_, R0, R1 ) * scalar;
            odres_  = repeat( mat_, R0, R1 ) * scalar;
            sres_   = repeat( mat_, R0, R1 ) * scalar;
            osres_  = repeat( mat_, R0, R1 ) * scalar;
            refres_ = repeat( refmat_, R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = repeat( omat_, R0, R1 ) * scalar;
            odres_  = repeat( omat_, R0, R1 ) * scalar;
            sres_   = repeat( omat_, R0, R1 ) * scalar;
            osres_  = repeat( omat_, R0, R1 ) * scalar;
            refres_ = repeat( refmat_, R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat operation with the given matrix (OP*s, compile time)
      {
         test_  = "Scaled repeat operation with the given matrix (OP*s, compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat<R0,R1>( mat_ ) * scalar;
            odres_  = repeat<R0,R1>( mat_ ) * scalar;
            sres_   = repeat<R0,R1>( mat_ ) * scalar;
            osres_  = repeat<R0,R1>( mat_ ) * scalar;
            refres_ = repeat<R0,R1>( refmat_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = repeat<R0,R1>( omat_ ) * scalar;
            odres_  = repeat<R0,R1>( omat_ ) * scalar;
            sres_   = repeat<R0,R1>( omat_ ) * scalar;
            osres_  = repeat<R0,R1>( omat_ ) * scalar;
            refres_ = repeat<R0,R1>( refmat_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat operation with evaluated matrix (OP*s, runtime)
      {
         test_  = "Scaled repeat operation with evaluated matrix (OP*s, runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat( eval( mat_ ), R0, R1 ) * scalar;
            odres_  = repeat( eval( mat_ ), R0, R1 ) * scalar;
            sres_   = repeat( eval( mat_ ), R0, R1 ) * scalar;
            osres_  = repeat( eval( mat_ ), R0, R1 ) * scalar;
            refres_ = repeat( eval( refmat_ ), R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = repeat( eval( omat_ ), R0, R1 ) * scalar;
            odres_  = repeat( eval( omat_ ), R0, R1 ) * scalar;
            sres_   = repeat( eval( omat_ ), R0, R1 ) * scalar;
            osres_  = repeat( eval( omat_ ), R0, R1 ) * scalar;
            refres_ = repeat( eval( refmat_ ), R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat operation with evaluated matrix (OP*s, compile time)
      {
         test_  = "Scaled repeat operation with the given matrix (OP*s, compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat<R0,R1>( eval( mat_ ) ) * scalar;
            odres_  = repeat<R0,R1>( eval( mat_ ) ) * scalar;
            sres_   = repeat<R0,R1>( eval( mat_ ) ) * scalar;
            osres_  = repeat<R0,R1>( eval( mat_ ) ) * scalar;
            refres_ = repeat<R0,R1>( eval( refmat_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = repeat<R0,R1>( eval( omat_ ) ) * scalar;
            odres_  = repeat<R0,R1>( eval( omat_ ) ) * scalar;
            sres_   = repeat<R0,R1>( eval( omat_ ) ) * scalar;
            osres_  = repeat<R0,R1>( eval( omat_ ) ) * scalar;
            refres_ = repeat<R0,R1>( eval( refmat_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Scaled repeat operation (OP/s)
      //=====================================================================================

      // Scaled repeat operation with the given matrix (OP/s, runtime)
      {
         test_  = "Scaled repeat operation with the given matrix (OP/s, runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat( mat_, R0, R1 ) / scalar;
            odres_  = repeat( mat_, R0, R1 ) / scalar;
            sres_   = repeat( mat_, R0, R1 ) / scalar;
            osres_  = repeat( mat_, R0, R1 ) / scalar;
            refres_ = repeat( refmat_, R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = repeat( omat_, R0, R1 ) / scalar;
            odres_  = repeat( omat_, R0, R1 ) / scalar;
            sres_   = repeat( omat_, R0, R1 ) / scalar;
            osres_  = repeat( omat_, R0, R1 ) / scalar;
            refres_ = repeat( refmat_, R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat operation with the given matrix (OP/s, compile time)
      {
         test_  = "Scaled repeat operation with the given matrix (OP/s, compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat<R0,R1>( mat_ ) / scalar;
            odres_  = repeat<R0,R1>( mat_ ) / scalar;
            sres_   = repeat<R0,R1>( mat_ ) / scalar;
            osres_  = repeat<R0,R1>( mat_ ) / scalar;
            refres_ = repeat<R0,R1>( refmat_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = repeat<R0,R1>( omat_ ) / scalar;
            odres_  = repeat<R0,R1>( omat_ ) / scalar;
            sres_   = repeat<R0,R1>( omat_ ) / scalar;
            osres_  = repeat<R0,R1>( omat_ ) / scalar;
            refres_ = repeat<R0,R1>( refmat_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat operation with evaluated matrix (OP/s, runtime)
      {
         test_  = "Scaled repeat operation with evaluated matrix (OP/s, runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat( eval( mat_ ), R0, R1 ) / scalar;
            odres_  = repeat( eval( mat_ ), R0, R1 ) / scalar;
            sres_   = repeat( eval( mat_ ), R0, R1 ) / scalar;
            osres_  = repeat( eval( mat_ ), R0, R1 ) / scalar;
            refres_ = repeat( eval( refmat_ ), R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = repeat( eval( omat_ ), R0, R1 ) / scalar;
            odres_  = repeat( eval( omat_ ), R0, R1 ) / scalar;
            sres_   = repeat( eval( omat_ ), R0, R1 ) / scalar;
            osres_  = repeat( eval( omat_ ), R0, R1 ) / scalar;
            refres_ = repeat( eval( refmat_ ), R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat operation with evaluated matrix (OP/s, compile time)
      {
         test_  = "Scaled repeat operation with the given matrix (OP/s, compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            dres_   = repeat<R0,R1>( eval( mat_ ) ) / scalar;
            odres_  = repeat<R0,R1>( eval( mat_ ) ) / scalar;
            sres_   = repeat<R0,R1>( eval( mat_ ) ) / scalar;
            osres_  = repeat<R0,R1>( eval( mat_ ) ) / scalar;
            refres_ = repeat<R0,R1>( eval( refmat_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   = repeat<R0,R1>( eval( omat_ ) ) / scalar;
            odres_  = repeat<R0,R1>( eval( omat_ ) ) / scalar;
            sres_   = repeat<R0,R1>( eval( omat_ ) ) / scalar;
            osres_  = repeat<R0,R1>( eval( omat_ ) ) / scalar;
            refres_ = repeat<R0,R1>( eval( refmat_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Scaled repeat with addition assignment (s*OP)
      //=====================================================================================

      // Scaled repeat with addition assignment with the given matrix (s*OP, runtime)
      {
         test_  = "Scaled repeat with addition assignment with the given matrix (s*OP, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += scalar * repeat( mat_, R0, R1 );
            odres_  += scalar * repeat( mat_, R0, R1 );
            sres_   += scalar * repeat( mat_, R0, R1 );
            osres_  += scalar * repeat( mat_, R0, R1 );
            refres_ += scalar * repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += scalar * repeat( omat_, R0, R1 );
            odres_  += scalar * repeat( omat_, R0, R1 );
            sres_   += scalar * repeat( omat_, R0, R1 );
            osres_  += scalar * repeat( omat_, R0, R1 );
            refres_ += scalar * repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with addition assignment with the given matrix (s*OP, compile time)
      {
         test_  = "Scaled repeat with addition assignment with the given matrix (s*OP, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += scalar * repeat<R0,R1>( mat_ );
            odres_  += scalar * repeat<R0,R1>( mat_ );
            sres_   += scalar * repeat<R0,R1>( mat_ );
            osres_  += scalar * repeat<R0,R1>( mat_ );
            refres_ += scalar * repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += scalar * repeat<R0,R1>( omat_ );
            odres_  += scalar * repeat<R0,R1>( omat_ );
            sres_   += scalar * repeat<R0,R1>( omat_ );
            osres_  += scalar * repeat<R0,R1>( omat_ );
            refres_ += scalar * repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with addition assignment with evaluated matrix (s*OP, runtime)
      {
         test_  = "Scaled repeat with addition assignment with evaluated matrix (s*OP, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += scalar * repeat( eval( mat_ ), R0, R1 );
            odres_  += scalar * repeat( eval( mat_ ), R0, R1 );
            sres_   += scalar * repeat( eval( mat_ ), R0, R1 );
            osres_  += scalar * repeat( eval( mat_ ), R0, R1 );
            refres_ += scalar * repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += scalar * repeat( eval( omat_ ), R0, R1 );
            odres_  += scalar * repeat( eval( omat_ ), R0, R1 );
            sres_   += scalar * repeat( eval( omat_ ), R0, R1 );
            osres_  += scalar * repeat( eval( omat_ ), R0, R1 );
            refres_ += scalar * repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with addition assignment with evaluated matrix (s*OP, compile time)
      {
         test_  = "Scaled repeat with addition assignment with the given matrix (s*OP, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += scalar * repeat<R0,R1>( eval( mat_ ) );
            odres_  += scalar * repeat<R0,R1>( eval( mat_ ) );
            sres_   += scalar * repeat<R0,R1>( eval( mat_ ) );
            osres_  += scalar * repeat<R0,R1>( eval( mat_ ) );
            refres_ += scalar * repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += scalar * repeat<R0,R1>( eval( omat_ ) );
            odres_  += scalar * repeat<R0,R1>( eval( omat_ ) );
            sres_   += scalar * repeat<R0,R1>( eval( omat_ ) );
            osres_  += scalar * repeat<R0,R1>( eval( omat_ ) );
            refres_ += scalar * repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Scaled repeat with addition assignment (OP*s)
      //=====================================================================================

      // Scaled repeat with addition assignment with the given matrix (OP*s, runtime)
      {
         test_  = "Scaled repeat with addition assignment with the given matrix (OP*s, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat( mat_, R0, R1 ) * scalar;
            odres_  += repeat( mat_, R0, R1 ) * scalar;
            sres_   += repeat( mat_, R0, R1 ) * scalar;
            osres_  += repeat( mat_, R0, R1 ) * scalar;
            refres_ += repeat( refmat_, R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += repeat( omat_, R0, R1 ) * scalar;
            odres_  += repeat( omat_, R0, R1 ) * scalar;
            sres_   += repeat( omat_, R0, R1 ) * scalar;
            osres_  += repeat( omat_, R0, R1 ) * scalar;
            refres_ += repeat( refmat_, R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with addition assignment with the given matrix (OP*s, compile time)
      {
         test_  = "Scaled repeat with addition assignment with the given matrix (OP*s, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat<R0,R1>( mat_ ) * scalar;
            odres_  += repeat<R0,R1>( mat_ ) * scalar;
            sres_   += repeat<R0,R1>( mat_ ) * scalar;
            osres_  += repeat<R0,R1>( mat_ ) * scalar;
            refres_ += repeat<R0,R1>( refmat_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += repeat<R0,R1>( omat_ ) * scalar;
            odres_  += repeat<R0,R1>( omat_ ) * scalar;
            sres_   += repeat<R0,R1>( omat_ ) * scalar;
            osres_  += repeat<R0,R1>( omat_ ) * scalar;
            refres_ += repeat<R0,R1>( refmat_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with addition assignment with evaluated matrix (OP*s, runtime)
      {
         test_  = "Scaled repeat with addition assignment with evaluated matrix (OP*s, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat( eval( mat_ ), R0, R1 ) * scalar;
            odres_  += repeat( eval( mat_ ), R0, R1 ) * scalar;
            sres_   += repeat( eval( mat_ ), R0, R1 ) * scalar;
            osres_  += repeat( eval( mat_ ), R0, R1 ) * scalar;
            refres_ += repeat( eval( refmat_ ), R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += repeat( eval( omat_ ), R0, R1 ) * scalar;
            odres_  += repeat( eval( omat_ ), R0, R1 ) * scalar;
            sres_   += repeat( eval( omat_ ), R0, R1 ) * scalar;
            osres_  += repeat( eval( omat_ ), R0, R1 ) * scalar;
            refres_ += repeat( eval( refmat_ ), R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with addition assignment with evaluated matrix (OP*s, compile time)
      {
         test_  = "Scaled repeat with addition assignment with the given matrix (OP*s, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat<R0,R1>( eval( mat_ ) ) * scalar;
            odres_  += repeat<R0,R1>( eval( mat_ ) ) * scalar;
            sres_   += repeat<R0,R1>( eval( mat_ ) ) * scalar;
            osres_  += repeat<R0,R1>( eval( mat_ ) ) * scalar;
            refres_ += repeat<R0,R1>( eval( refmat_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += repeat<R0,R1>( eval( omat_ ) ) * scalar;
            odres_  += repeat<R0,R1>( eval( omat_ ) ) * scalar;
            sres_   += repeat<R0,R1>( eval( omat_ ) ) * scalar;
            osres_  += repeat<R0,R1>( eval( omat_ ) ) * scalar;
            refres_ += repeat<R0,R1>( eval( refmat_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Scaled repeat with addition assignment (OP/s)
      //=====================================================================================

      // Scaled repeat with addition assignment with the given matrix (OP/s, runtime)
      {
         test_  = "Scaled repeat with addition assignment with the given matrix (OP/s, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat( mat_, R0, R1 ) / scalar;
            odres_  += repeat( mat_, R0, R1 ) / scalar;
            sres_   += repeat( mat_, R0, R1 ) / scalar;
            osres_  += repeat( mat_, R0, R1 ) / scalar;
            refres_ += repeat( refmat_, R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += repeat( omat_, R0, R1 ) / scalar;
            odres_  += repeat( omat_, R0, R1 ) / scalar;
            sres_   += repeat( omat_, R0, R1 ) / scalar;
            osres_  += repeat( omat_, R0, R1 ) / scalar;
            refres_ += repeat( refmat_, R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with addition assignment with the given matrix (OP/s, compile time)
      {
         test_  = "Scaled repeat with addition assignment with the given matrix (OP/s, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat<R0,R1>( mat_ ) / scalar;
            odres_  += repeat<R0,R1>( mat_ ) / scalar;
            sres_   += repeat<R0,R1>( mat_ ) / scalar;
            osres_  += repeat<R0,R1>( mat_ ) / scalar;
            refres_ += repeat<R0,R1>( refmat_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += repeat<R0,R1>( omat_ ) / scalar;
            odres_  += repeat<R0,R1>( omat_ ) / scalar;
            sres_   += repeat<R0,R1>( omat_ ) / scalar;
            osres_  += repeat<R0,R1>( omat_ ) / scalar;
            refres_ += repeat<R0,R1>( refmat_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with addition assignment with evaluated matrix (OP/s, runtime)
      {
         test_  = "Scaled repeat with addition assignment with evaluated matrix (OP/s, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat( eval( mat_ ), R0, R1 ) / scalar;
            odres_  += repeat( eval( mat_ ), R0, R1 ) / scalar;
            sres_   += repeat( eval( mat_ ), R0, R1 ) / scalar;
            osres_  += repeat( eval( mat_ ), R0, R1 ) / scalar;
            refres_ += repeat( eval( refmat_ ), R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += repeat( eval( omat_ ), R0, R1 ) / scalar;
            odres_  += repeat( eval( omat_ ), R0, R1 ) / scalar;
            sres_   += repeat( eval( omat_ ), R0, R1 ) / scalar;
            osres_  += repeat( eval( omat_ ), R0, R1 ) / scalar;
            refres_ += repeat( eval( refmat_ ), R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with addition assignment with evaluated matrix (OP/s, compile time)
      {
         test_  = "Scaled repeat with addition assignment with the given matrix (OP/s, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += repeat<R0,R1>( eval( mat_ ) ) / scalar;
            odres_  += repeat<R0,R1>( eval( mat_ ) ) / scalar;
            sres_   += repeat<R0,R1>( eval( mat_ ) ) / scalar;
            osres_  += repeat<R0,R1>( eval( mat_ ) ) / scalar;
            refres_ += repeat<R0,R1>( eval( refmat_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   += repeat<R0,R1>( eval( omat_ ) ) / scalar;
            odres_  += repeat<R0,R1>( eval( omat_ ) ) / scalar;
            sres_   += repeat<R0,R1>( eval( omat_ ) ) / scalar;
            osres_  += repeat<R0,R1>( eval( omat_ ) ) / scalar;
            refres_ += repeat<R0,R1>( eval( refmat_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Scaled repeat with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled repeat with subtraction assignment with the given matrix (s*OP, runtime)
      {
         test_  = "Scaled repeat with subtraction assignment with the given matrix (s*OP, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= scalar * repeat( mat_, R0, R1 );
            odres_  -= scalar * repeat( mat_, R0, R1 );
            sres_   -= scalar * repeat( mat_, R0, R1 );
            osres_  -= scalar * repeat( mat_, R0, R1 );
            refres_ -= scalar * repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= scalar * repeat( omat_, R0, R1 );
            odres_  -= scalar * repeat( omat_, R0, R1 );
            sres_   -= scalar * repeat( omat_, R0, R1 );
            osres_  -= scalar * repeat( omat_, R0, R1 );
            refres_ -= scalar * repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with subtraction assignment with the given matrix (s*OP, compile time)
      {
         test_  = "Scaled repeat with subtraction assignment with the given matrix (s*OP, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= scalar * repeat<R0,R1>( mat_ );
            odres_  -= scalar * repeat<R0,R1>( mat_ );
            sres_   -= scalar * repeat<R0,R1>( mat_ );
            osres_  -= scalar * repeat<R0,R1>( mat_ );
            refres_ -= scalar * repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= scalar * repeat<R0,R1>( omat_ );
            odres_  -= scalar * repeat<R0,R1>( omat_ );
            sres_   -= scalar * repeat<R0,R1>( omat_ );
            osres_  -= scalar * repeat<R0,R1>( omat_ );
            refres_ -= scalar * repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with subtraction assignment with evaluated matrix (s*OP, runtime)
      {
         test_  = "Scaled repeat with subtraction assignment with evaluated matrix (s*OP, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= scalar * repeat( eval( mat_ ), R0, R1 );
            odres_  -= scalar * repeat( eval( mat_ ), R0, R1 );
            sres_   -= scalar * repeat( eval( mat_ ), R0, R1 );
            osres_  -= scalar * repeat( eval( mat_ ), R0, R1 );
            refres_ -= scalar * repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= scalar * repeat( eval( omat_ ), R0, R1 );
            odres_  -= scalar * repeat( eval( omat_ ), R0, R1 );
            sres_   -= scalar * repeat( eval( omat_ ), R0, R1 );
            osres_  -= scalar * repeat( eval( omat_ ), R0, R1 );
            refres_ -= scalar * repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with subtraction assignment with evaluated matrix (s*OP, compile time)
      {
         test_  = "Scaled repeat with subtraction assignment with the given matrix (s*OP, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= scalar * repeat<R0,R1>( eval( mat_ ) );
            odres_  -= scalar * repeat<R0,R1>( eval( mat_ ) );
            sres_   -= scalar * repeat<R0,R1>( eval( mat_ ) );
            osres_  -= scalar * repeat<R0,R1>( eval( mat_ ) );
            refres_ -= scalar * repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= scalar * repeat<R0,R1>( eval( omat_ ) );
            odres_  -= scalar * repeat<R0,R1>( eval( omat_ ) );
            sres_   -= scalar * repeat<R0,R1>( eval( omat_ ) );
            osres_  -= scalar * repeat<R0,R1>( eval( omat_ ) );
            refres_ -= scalar * repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Scaled repeat with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled repeat with subtraction assignment with the given matrix (OP*s, runtime)
      {
         test_  = "Scaled repeat with subtraction assignment with the given matrix (OP*s, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat( mat_, R0, R1 ) * scalar;
            odres_  -= repeat( mat_, R0, R1 ) * scalar;
            sres_   -= repeat( mat_, R0, R1 ) * scalar;
            osres_  -= repeat( mat_, R0, R1 ) * scalar;
            refres_ -= repeat( refmat_, R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= repeat( omat_, R0, R1 ) * scalar;
            odres_  -= repeat( omat_, R0, R1 ) * scalar;
            sres_   -= repeat( omat_, R0, R1 ) * scalar;
            osres_  -= repeat( omat_, R0, R1 ) * scalar;
            refres_ -= repeat( refmat_, R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with subtraction assignment with the given matrix (OP*s, compile time)
      {
         test_  = "Scaled repeat with subtraction assignment with the given matrix (OP*s, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat<R0,R1>( mat_ ) * scalar;
            odres_  -= repeat<R0,R1>( mat_ ) * scalar;
            sres_   -= repeat<R0,R1>( mat_ ) * scalar;
            osres_  -= repeat<R0,R1>( mat_ ) * scalar;
            refres_ -= repeat<R0,R1>( refmat_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= repeat<R0,R1>( omat_ ) * scalar;
            odres_  -= repeat<R0,R1>( omat_ ) * scalar;
            sres_   -= repeat<R0,R1>( omat_ ) * scalar;
            osres_  -= repeat<R0,R1>( omat_ ) * scalar;
            refres_ -= repeat<R0,R1>( refmat_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with subtraction assignment with evaluated matrix (OP*s, runtime)
      {
         test_  = "Scaled repeat with subtraction assignment with evaluated matrix (OP*s, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat( eval( mat_ ), R0, R1 ) * scalar;
            odres_  -= repeat( eval( mat_ ), R0, R1 ) * scalar;
            sres_   -= repeat( eval( mat_ ), R0, R1 ) * scalar;
            osres_  -= repeat( eval( mat_ ), R0, R1 ) * scalar;
            refres_ -= repeat( eval( refmat_ ), R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= repeat( eval( omat_ ), R0, R1 ) * scalar;
            odres_  -= repeat( eval( omat_ ), R0, R1 ) * scalar;
            sres_   -= repeat( eval( omat_ ), R0, R1 ) * scalar;
            osres_  -= repeat( eval( omat_ ), R0, R1 ) * scalar;
            refres_ -= repeat( eval( refmat_ ), R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with subtraction assignment with evaluated matrix (OP*s, compile time)
      {
         test_  = "Scaled repeat with subtraction assignment with the given matrix (OP*s, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat<R0,R1>( eval( mat_ ) ) * scalar;
            odres_  -= repeat<R0,R1>( eval( mat_ ) ) * scalar;
            sres_   -= repeat<R0,R1>( eval( mat_ ) ) * scalar;
            osres_  -= repeat<R0,R1>( eval( mat_ ) ) * scalar;
            refres_ -= repeat<R0,R1>( eval( refmat_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= repeat<R0,R1>( eval( omat_ ) ) * scalar;
            odres_  -= repeat<R0,R1>( eval( omat_ ) ) * scalar;
            sres_   -= repeat<R0,R1>( eval( omat_ ) ) * scalar;
            osres_  -= repeat<R0,R1>( eval( omat_ ) ) * scalar;
            refres_ -= repeat<R0,R1>( eval( refmat_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Scaled repeat with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled repeat with subtraction assignment with the given matrix (OP/s, runtime)
      {
         test_  = "Scaled repeat with subtraction assignment with the given matrix (OP/s, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat( mat_, R0, R1 ) / scalar;
            odres_  -= repeat( mat_, R0, R1 ) / scalar;
            sres_   -= repeat( mat_, R0, R1 ) / scalar;
            osres_  -= repeat( mat_, R0, R1 ) / scalar;
            refres_ -= repeat( refmat_, R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= repeat( omat_, R0, R1 ) / scalar;
            odres_  -= repeat( omat_, R0, R1 ) / scalar;
            sres_   -= repeat( omat_, R0, R1 ) / scalar;
            osres_  -= repeat( omat_, R0, R1 ) / scalar;
            refres_ -= repeat( refmat_, R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with subtraction assignment with the given matrix (OP/s, compile time)
      {
         test_  = "Scaled repeat with subtraction assignment with the given matrix (OP/s, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat<R0,R1>( mat_ ) / scalar;
            odres_  -= repeat<R0,R1>( mat_ ) / scalar;
            sres_   -= repeat<R0,R1>( mat_ ) / scalar;
            osres_  -= repeat<R0,R1>( mat_ ) / scalar;
            refres_ -= repeat<R0,R1>( refmat_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= repeat<R0,R1>( omat_ ) / scalar;
            odres_  -= repeat<R0,R1>( omat_ ) / scalar;
            sres_   -= repeat<R0,R1>( omat_ ) / scalar;
            osres_  -= repeat<R0,R1>( omat_ ) / scalar;
            refres_ -= repeat<R0,R1>( refmat_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with subtraction assignment with evaluated matrix (OP/s, runtime)
      {
         test_  = "Scaled repeat with subtraction assignment with evaluated matrix (OP/s, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat( eval( mat_ ), R0, R1 ) / scalar;
            odres_  -= repeat( eval( mat_ ), R0, R1 ) / scalar;
            sres_   -= repeat( eval( mat_ ), R0, R1 ) / scalar;
            osres_  -= repeat( eval( mat_ ), R0, R1 ) / scalar;
            refres_ -= repeat( eval( refmat_ ), R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= repeat( eval( omat_ ), R0, R1 ) / scalar;
            odres_  -= repeat( eval( omat_ ), R0, R1 ) / scalar;
            sres_   -= repeat( eval( omat_ ), R0, R1 ) / scalar;
            osres_  -= repeat( eval( omat_ ), R0, R1 ) / scalar;
            refres_ -= repeat( eval( refmat_ ), R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with subtraction assignment with evaluated matrix (OP/s, compile time)
      {
         test_  = "Scaled repeat with subtraction assignment with the given matrix (OP/s, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= repeat<R0,R1>( eval( mat_ ) ) / scalar;
            odres_  -= repeat<R0,R1>( eval( mat_ ) ) / scalar;
            sres_   -= repeat<R0,R1>( eval( mat_ ) ) / scalar;
            osres_  -= repeat<R0,R1>( eval( mat_ ) ) / scalar;
            refres_ -= repeat<R0,R1>( eval( refmat_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   -= repeat<R0,R1>( eval( omat_ ) ) / scalar;
            odres_  -= repeat<R0,R1>( eval( omat_ ) ) / scalar;
            sres_   -= repeat<R0,R1>( eval( omat_ ) ) / scalar;
            osres_  -= repeat<R0,R1>( eval( omat_ ) ) / scalar;
            refres_ -= repeat<R0,R1>( eval( refmat_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Scaled repeat with Schur product assignment (s*OP)
      //=====================================================================================

      // Scaled repeat with Schur product assignment with the given matrix (s*OP, runtime)
      {
         test_  = "Scaled repeat with Schur product assignment with the given matrix (s*OP, runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= scalar * repeat( mat_, R0, R1 );
            odres_  %= scalar * repeat( mat_, R0, R1 );
            sres_   %= scalar * repeat( mat_, R0, R1 );
            osres_  %= scalar * repeat( mat_, R0, R1 );
            refres_ %= scalar * repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= scalar * repeat( omat_, R0, R1 );
            odres_  %= scalar * repeat( omat_, R0, R1 );
            sres_   %= scalar * repeat( omat_, R0, R1 );
            osres_  %= scalar * repeat( omat_, R0, R1 );
            refres_ %= scalar * repeat( refmat_, R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with Schur product assignment with the given matrix (s*OP, compile time)
      {
         test_  = "Scaled repeat with Schur product assignment with the given matrix (s*OP, compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= scalar * repeat<R0,R1>( mat_ );
            odres_  %= scalar * repeat<R0,R1>( mat_ );
            sres_   %= scalar * repeat<R0,R1>( mat_ );
            osres_  %= scalar * repeat<R0,R1>( mat_ );
            refres_ %= scalar * repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= scalar * repeat<R0,R1>( omat_ );
            odres_  %= scalar * repeat<R0,R1>( omat_ );
            sres_   %= scalar * repeat<R0,R1>( omat_ );
            osres_  %= scalar * repeat<R0,R1>( omat_ );
            refres_ %= scalar * repeat<R0,R1>( refmat_ );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with Schur product assignment with evaluated matrix (s*OP, runtime)
      {
         test_  = "Scaled repeat with Schur product assignment with evaluated matrix (s*OP, runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= scalar * repeat( eval( mat_ ), R0, R1 );
            odres_  %= scalar * repeat( eval( mat_ ), R0, R1 );
            sres_   %= scalar * repeat( eval( mat_ ), R0, R1 );
            osres_  %= scalar * repeat( eval( mat_ ), R0, R1 );
            refres_ %= scalar * repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= scalar * repeat( eval( omat_ ), R0, R1 );
            odres_  %= scalar * repeat( eval( omat_ ), R0, R1 );
            sres_   %= scalar * repeat( eval( omat_ ), R0, R1 );
            osres_  %= scalar * repeat( eval( omat_ ), R0, R1 );
            refres_ %= scalar * repeat( eval( refmat_ ), R0, R1 );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with Schur product assignment with evaluated matrix (s*OP, compile time)
      {
         test_  = "Scaled repeat with Schur product assignment with the given matrix (s*OP, compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= scalar * repeat<R0,R1>( eval( mat_ ) );
            odres_  %= scalar * repeat<R0,R1>( eval( mat_ ) );
            sres_   %= scalar * repeat<R0,R1>( eval( mat_ ) );
            osres_  %= scalar * repeat<R0,R1>( eval( mat_ ) );
            refres_ %= scalar * repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= scalar * repeat<R0,R1>( eval( omat_ ) );
            odres_  %= scalar * repeat<R0,R1>( eval( omat_ ) );
            sres_   %= scalar * repeat<R0,R1>( eval( omat_ ) );
            osres_  %= scalar * repeat<R0,R1>( eval( omat_ ) );
            refres_ %= scalar * repeat<R0,R1>( eval( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Scaled repeat with Schur product assignment (OP*s)
      //=====================================================================================

      // Scaled repeat with Schur product assignment with the given matrix (OP*s, runtime)
      {
         test_  = "Scaled repeat with Schur product assignment with the given matrix (OP*s, runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= repeat( mat_, R0, R1 ) * scalar;
            odres_  %= repeat( mat_, R0, R1 ) * scalar;
            sres_   %= repeat( mat_, R0, R1 ) * scalar;
            osres_  %= repeat( mat_, R0, R1 ) * scalar;
            refres_ %= repeat( refmat_, R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= repeat( omat_, R0, R1 ) * scalar;
            odres_  %= repeat( omat_, R0, R1 ) * scalar;
            sres_   %= repeat( omat_, R0, R1 ) * scalar;
            osres_  %= repeat( omat_, R0, R1 ) * scalar;
            refres_ %= repeat( refmat_, R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with Schur product assignment with the given matrix (OP*s, compile time)
      {
         test_  = "Scaled repeat with Schur product assignment with the given matrix (OP*s, compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= repeat<R0,R1>( mat_ ) * scalar;
            odres_  %= repeat<R0,R1>( mat_ ) * scalar;
            sres_   %= repeat<R0,R1>( mat_ ) * scalar;
            osres_  %= repeat<R0,R1>( mat_ ) * scalar;
            refres_ %= repeat<R0,R1>( refmat_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= repeat<R0,R1>( omat_ ) * scalar;
            odres_  %= repeat<R0,R1>( omat_ ) * scalar;
            sres_   %= repeat<R0,R1>( omat_ ) * scalar;
            osres_  %= repeat<R0,R1>( omat_ ) * scalar;
            refres_ %= repeat<R0,R1>( refmat_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with Schur product assignment with evaluated matrix (OP*s, runtime)
      {
         test_  = "Scaled repeat with Schur product assignment with evaluated matrix (OP*s, runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= repeat( eval( mat_ ), R0, R1 ) * scalar;
            odres_  %= repeat( eval( mat_ ), R0, R1 ) * scalar;
            sres_   %= repeat( eval( mat_ ), R0, R1 ) * scalar;
            osres_  %= repeat( eval( mat_ ), R0, R1 ) * scalar;
            refres_ %= repeat( eval( refmat_ ), R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= repeat( eval( omat_ ), R0, R1 ) * scalar;
            odres_  %= repeat( eval( omat_ ), R0, R1 ) * scalar;
            sres_   %= repeat( eval( omat_ ), R0, R1 ) * scalar;
            osres_  %= repeat( eval( omat_ ), R0, R1 ) * scalar;
            refres_ %= repeat( eval( refmat_ ), R0, R1 ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with Schur product assignment with evaluated matrix (OP*s, compile time)
      {
         test_  = "Scaled repeat with Schur product assignment with the given matrix (OP*s, compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= repeat<R0,R1>( eval( mat_ ) ) * scalar;
            odres_  %= repeat<R0,R1>( eval( mat_ ) ) * scalar;
            sres_   %= repeat<R0,R1>( eval( mat_ ) ) * scalar;
            osres_  %= repeat<R0,R1>( eval( mat_ ) ) * scalar;
            refres_ %= repeat<R0,R1>( eval( refmat_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= repeat<R0,R1>( eval( omat_ ) ) * scalar;
            odres_  %= repeat<R0,R1>( eval( omat_ ) ) * scalar;
            sres_   %= repeat<R0,R1>( eval( omat_ ) ) * scalar;
            osres_  %= repeat<R0,R1>( eval( omat_ ) ) * scalar;
            refres_ %= repeat<R0,R1>( eval( refmat_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Scaled repeat with Schur product assignment (OP/s)
      //=====================================================================================

      // Scaled repeat with Schur product assignment with the given matrix (OP/s, runtime)
      {
         test_  = "Scaled repeat with Schur product assignment with the given matrix (OP/s, runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= repeat( mat_, R0, R1 ) / scalar;
            odres_  %= repeat( mat_, R0, R1 ) / scalar;
            sres_   %= repeat( mat_, R0, R1 ) / scalar;
            osres_  %= repeat( mat_, R0, R1 ) / scalar;
            refres_ %= repeat( refmat_, R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= repeat( omat_, R0, R1 ) / scalar;
            odres_  %= repeat( omat_, R0, R1 ) / scalar;
            sres_   %= repeat( omat_, R0, R1 ) / scalar;
            osres_  %= repeat( omat_, R0, R1 ) / scalar;
            refres_ %= repeat( refmat_, R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with Schur product assignment with the given matrix (OP/s, compile time)
      {
         test_  = "Scaled repeat with Schur product assignment with the given matrix (OP/s, compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= repeat<R0,R1>( mat_ ) / scalar;
            odres_  %= repeat<R0,R1>( mat_ ) / scalar;
            sres_   %= repeat<R0,R1>( mat_ ) / scalar;
            osres_  %= repeat<R0,R1>( mat_ ) / scalar;
            refres_ %= repeat<R0,R1>( refmat_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= repeat<R0,R1>( omat_ ) / scalar;
            odres_  %= repeat<R0,R1>( omat_ ) / scalar;
            sres_   %= repeat<R0,R1>( omat_ ) / scalar;
            osres_  %= repeat<R0,R1>( omat_ ) / scalar;
            refres_ %= repeat<R0,R1>( refmat_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with Schur product assignment with evaluated matrix (OP/s, runtime)
      {
         test_  = "Scaled repeat with Schur product assignment with evaluated matrix (OP/s, runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= repeat( eval( mat_ ), R0, R1 ) / scalar;
            odres_  %= repeat( eval( mat_ ), R0, R1 ) / scalar;
            sres_   %= repeat( eval( mat_ ), R0, R1 ) / scalar;
            osres_  %= repeat( eval( mat_ ), R0, R1 ) / scalar;
            refres_ %= repeat( eval( refmat_ ), R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= repeat( eval( omat_ ), R0, R1 ) / scalar;
            odres_  %= repeat( eval( omat_ ), R0, R1 ) / scalar;
            sres_   %= repeat( eval( omat_ ), R0, R1 ) / scalar;
            osres_  %= repeat( eval( omat_ ), R0, R1 ) / scalar;
            refres_ %= repeat( eval( refmat_ ), R0, R1 ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Scaled repeat with Schur product assignment with evaluated matrix (OP/s, compile time)
      {
         test_  = "Scaled repeat with Schur product assignment with the given matrix (OP/s, compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= repeat<R0,R1>( eval( mat_ ) ) / scalar;
            odres_  %= repeat<R0,R1>( eval( mat_ ) ) / scalar;
            sres_   %= repeat<R0,R1>( eval( mat_ ) ) / scalar;
            osres_  %= repeat<R0,R1>( eval( mat_ ) ) / scalar;
            refres_ %= repeat<R0,R1>( eval( refmat_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            dres_   %= repeat<R0,R1>( eval( omat_ ) ) / scalar;
            odres_  %= repeat<R0,R1>( eval( omat_ ) ) / scalar;
            sres_   %= repeat<R0,R1>( eval( omat_ ) ) / scalar;
            osres_  %= repeat<R0,R1>( eval( omat_ ) ) / scalar;
            refres_ %= repeat<R0,R1>( eval( refmat_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the transpose sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the transpose matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANS_OPERATION > 1 )
   {
      using blaze::repeat;


      //=====================================================================================
      // Transpose repeat operation
      //=====================================================================================

      // Transpose repeat operation with the given matrix (runtime)
      {
         test_  = "Transpose repeat operation with the given matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_  = trans( repeat( mat_, R0, R1 ) );
            todres_ = trans( repeat( mat_, R0, R1 ) );
            tsres_  = trans( repeat( mat_, R0, R1 ) );
            tosres_ = trans( repeat( mat_, R0, R1 ) );
            refres_ = trans( repeat( refmat_, R0, R1 ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_  = trans( repeat( omat_, R0, R1 ) );
            todres_ = trans( repeat( omat_, R0, R1 ) );
            tsres_  = trans( repeat( omat_, R0, R1 ) );
            tosres_ = trans( repeat( omat_, R0, R1 ) );
            refres_ = trans( repeat( refmat_, R0, R1 ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkTransposeResults<OMT>();
      }

      // Transpose repeat operation with the given matrix (compile time)
      {
         test_  = "Transpose repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_  = trans( repeat<R0,R1>( mat_ ) );
            todres_ = trans( repeat<R0,R1>( mat_ ) );
            tsres_  = trans( repeat<R0,R1>( mat_ ) );
            tosres_ = trans( repeat<R0,R1>( mat_ ) );
            refres_ = trans( repeat<R0,R1>( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_  = trans( repeat<R0,R1>( omat_ ) );
            todres_ = trans( repeat<R0,R1>( omat_ ) );
            tsres_  = trans( repeat<R0,R1>( omat_ ) );
            tosres_ = trans( repeat<R0,R1>( omat_ ) );
            refres_ = trans( repeat<R0,R1>( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkTransposeResults<OMT>();
      }

      // Transpose repeat operation with evaluated matrix (runtime)
      {
         test_  = "Transpose repeat operation with evaluated matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_  = trans( repeat( eval( mat_ ), R0, R1 ) );
            todres_ = trans( repeat( eval( mat_ ), R0, R1 ) );
            tsres_  = trans( repeat( eval( mat_ ), R0, R1 ) );
            tosres_ = trans( repeat( eval( mat_ ), R0, R1 ) );
            refres_ = trans( repeat( eval( refmat_ ), R0, R1 ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_  = trans( repeat( eval( omat_ ), R0, R1 ) );
            todres_ = trans( repeat( eval( omat_ ), R0, R1 ) );
            tsres_  = trans( repeat( eval( omat_ ), R0, R1 ) );
            tosres_ = trans( repeat( eval( omat_ ), R0, R1 ) );
            refres_ = trans( repeat( eval( refmat_ ), R0, R1 ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkTransposeResults<OMT>();
      }

      // Transpose repeat operation with evaluated matrix (compile time)
      {
         test_  = "Transpose repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_  = trans( repeat<R0,R1>( eval( mat_ ) ) );
            todres_ = trans( repeat<R0,R1>( eval( mat_ ) ) );
            tsres_  = trans( repeat<R0,R1>( eval( mat_ ) ) );
            tosres_ = trans( repeat<R0,R1>( eval( mat_ ) ) );
            refres_ = trans( repeat<R0,R1>( eval( refmat_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_  = trans( repeat<R0,R1>( eval( omat_ ) ) );
            todres_ = trans( repeat<R0,R1>( eval( omat_ ) ) );
            tsres_  = trans( repeat<R0,R1>( eval( omat_ ) ) );
            tosres_ = trans( repeat<R0,R1>( eval( omat_ ) ) );
            refres_ = trans( repeat<R0,R1>( eval( refmat_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkTransposeResults<OMT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the conjugate transpose sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the conjugate transpose matrix repeat operation with plain assignment,
// addition assignment, subtraction assignment, multiplication assignment, and division assignment.
// In case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testCTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION > 1 )
   {
      using blaze::repeat;


      //=====================================================================================
      // Conjugate transpose repeat operation
      //=====================================================================================

      // Conjugate transpose repeat operation with the given matrix (runtime)
      {
         test_  = "Conjugate transpose repeat operation with the given matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_  = ctrans( repeat( mat_, R0, R1 ) );
            todres_ = ctrans( repeat( mat_, R0, R1 ) );
            tsres_  = ctrans( repeat( mat_, R0, R1 ) );
            tosres_ = ctrans( repeat( mat_, R0, R1 ) );
            refres_ = ctrans( repeat( refmat_, R0, R1 ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_  = ctrans( repeat( omat_, R0, R1 ) );
            todres_ = ctrans( repeat( omat_, R0, R1 ) );
            tsres_  = ctrans( repeat( omat_, R0, R1 ) );
            tosres_ = ctrans( repeat( omat_, R0, R1 ) );
            refres_ = ctrans( repeat( refmat_, R0, R1 ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkTransposeResults<OMT>();
      }

      // Conjugate transpose repeat operation with the given matrix (compile time)
      {
         test_  = "Conjugate transpose repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_  = ctrans( repeat<R0,R1>( mat_ ) );
            todres_ = ctrans( repeat<R0,R1>( mat_ ) );
            tsres_  = ctrans( repeat<R0,R1>( mat_ ) );
            tosres_ = ctrans( repeat<R0,R1>( mat_ ) );
            refres_ = ctrans( repeat<R0,R1>( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_  = ctrans( repeat<R0,R1>( omat_ ) );
            todres_ = ctrans( repeat<R0,R1>( omat_ ) );
            tsres_  = ctrans( repeat<R0,R1>( omat_ ) );
            tosres_ = ctrans( repeat<R0,R1>( omat_ ) );
            refres_ = ctrans( repeat<R0,R1>( refmat_ ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkTransposeResults<OMT>();
      }

      // Conjugate transpose repeat operation with evaluated matrix (runtime)
      {
         test_  = "Conjugate transpose repeat operation with evaluated matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_  = ctrans( repeat( eval( mat_ ), R0, R1 ) );
            todres_ = ctrans( repeat( eval( mat_ ), R0, R1 ) );
            tsres_  = ctrans( repeat( eval( mat_ ), R0, R1 ) );
            tosres_ = ctrans( repeat( eval( mat_ ), R0, R1 ) );
            refres_ = ctrans( repeat( eval( refmat_ ), R0, R1 ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_  = ctrans( repeat( eval( omat_ ), R0, R1 ) );
            todres_ = ctrans( repeat( eval( omat_ ), R0, R1 ) );
            tsres_  = ctrans( repeat( eval( omat_ ), R0, R1 ) );
            tosres_ = ctrans( repeat( eval( omat_ ), R0, R1 ) );
            refres_ = ctrans( repeat( eval( refmat_ ), R0, R1 ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkTransposeResults<OMT>();
      }

      // Conjugate transpose repeat operation with evaluated matrix (compile time)
      {
         test_  = "Conjugate transpose repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initTransposeResults();
            tdres_  = ctrans( repeat<R0,R1>( eval( mat_ ) ) );
            todres_ = ctrans( repeat<R0,R1>( eval( mat_ ) ) );
            tsres_  = ctrans( repeat<R0,R1>( eval( mat_ ) ) );
            tosres_ = ctrans( repeat<R0,R1>( eval( mat_ ) ) );
            refres_ = ctrans( repeat<R0,R1>( eval( refmat_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkTransposeResults<MT>();

         try {
            initTransposeResults();
            tdres_  = ctrans( repeat<R0,R1>( eval( omat_ ) ) );
            todres_ = ctrans( repeat<R0,R1>( eval( omat_ ) ) );
            tsres_  = ctrans( repeat<R0,R1>( eval( omat_ ) ) );
            tosres_ = ctrans( repeat<R0,R1>( eval( omat_ ) ) );
            refres_ = ctrans( repeat<R0,R1>( eval( refmat_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkTransposeResults<OMT>();
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the abs sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the abs matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testAbsOperation()
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
/*!\brief Testing the conjugate sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the conjugate matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testConjOperation()
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
/*!\brief Testing the \a real sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the \a real matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testRealOperation()
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
/*!\brief Testing the \a imag sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the \a imag matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testImagOperation()
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
/*!\brief Testing the evaluated sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the evaluated matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testEvalOperation()
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
/*!\brief Testing the serialized sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the serialized matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testSerialOperation()
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
/*!\brief Testing the non-aliased sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the non-aliased matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testNoAliasOperation()
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
/*!\brief Testing the non-SIMD sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the non-SIMD matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testNoSIMDOperation()
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
/*!\brief Testing the submatrix-wise sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the submatrix-wise matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the repeat operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testSubmatrixOperation()
{
#if BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION > 1 )
   {
      using blaze::repeat;


      if( mat_.rows() == 0UL || mat_.columns() == 0UL )
         return;


      //=====================================================================================
      // Submatrix-wise repeat operation
      //=====================================================================================

      // Submatrix-wise repeat operation with the given matrix (runtime)
      {
         test_  = "Submatrix-wise repeat operation with the given matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( repeat( refmat_, R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( repeat( refmat_, R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Submatrix-wise repeat operation with the given matrix (compile time)
      {
         test_  = "Submatrix-wise repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( repeat<R0,R1>( refmat_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( repeat<R0,R1>( refmat_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Submatrix-wise repeat operation with evaluated matrix (runtime)
      {
         test_  = "Submatrix-wise repeat operation with evaluated matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( repeat( eval( refmat_ ), R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( repeat( eval( refmat_ ), R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Submatrix-wise repeat operation with evaluated matrix (compile time)
      {
         test_  = "Submatrix-wise repeat operation with evaluated matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( repeat<R0,R1>( eval( refmat_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( repeat<R0,R1>( eval( refmat_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Submatrix-wise repeat with addition assignment
      //=====================================================================================

      // Submatrix-wise repeat with addition assignment with the given matrix (runtime)
      {
         test_  = "Submatrix-wise repeat with addition assignment with the given matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( repeat( refmat_, R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( repeat( refmat_, R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Submatrix-wise repeat with addition assignment with the given matrix (compile time)
      {
         test_  = "Submatrix-wise repeat with addition assignment with the given matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( repeat<R0,R1>( refmat_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( repeat<R0,R1>( refmat_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Submatrix-wise repeat with addition assignment with evaluated matrix (runtime)
      {
         test_  = "Submatrix-wise repeat with addition assignment with evaluated matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( repeat( eval( refmat_ ), R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( repeat( eval( refmat_ ), R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Submatrix-wise repeat with addition assignment with evaluated matrix (compile time)
      {
         test_  = "Submatrix-wise repeat with addition assignment with evaluated matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( repeat<R0,R1>( eval( refmat_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( repeat<R0,R1>( eval( refmat_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Submatrix-wise repeat with subtraction assignment
      //=====================================================================================

      // Submatrix-wise repeat with subtraction assignment with the given matrix (runtime)
      {
         test_  = "Submatrix-wise repeat with subtraction assignment with the given matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( repeat( refmat_, R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( repeat( refmat_, R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Submatrix-wise repeat with subtraction assignment with the given matrix (compile time)
      {
         test_  = "Submatrix-wise repeat with subtraction assignment with the given matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( repeat<R0,R1>( refmat_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( repeat<R0,R1>( refmat_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Submatrix-wise repeat with subtraction assignment with evaluated matrix (runtime)
      {
         test_  = "Submatrix-wise repeat with subtraction assignment with evaluated matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( repeat( eval( refmat_ ), R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( repeat( eval( refmat_ ), R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Submatrix-wise repeat with subtraction assignment with evaluated matrix (compile time)
      {
         test_  = "Submatrix-wise repeat with subtraction assignment with evaluated matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( repeat<R0,R1>( eval( refmat_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( repeat<R0,R1>( eval( refmat_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Submatrix-wise repeat with Schur product assignment
      //=====================================================================================

      // Submatrix-wise repeat with Schur product assignment with the given matrix (runtime)
      {
         test_  = "Submatrix-wise repeat with Schur product assignment with the given matrix (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( repeat( mat_   , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( repeat( refmat_, R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( repeat( omat_  , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( repeat( refmat_, R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Submatrix-wise repeat with Schur product assignment with the given matrix (compile time)
      {
         test_  = "Submatrix-wise repeat with Schur product assignment with the given matrix (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( repeat<R0,R1>( mat_    ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( repeat<R0,R1>( refmat_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( repeat<R0,R1>( omat_   ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( repeat<R0,R1>( refmat_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Submatrix-wise repeat with Schur product assignment with evaluated matrix (runtime)
      {
         test_  = "Submatrix-wise repeat with Schur product assignment with evaluated matrix (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( repeat( eval( mat_ )   , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( repeat( eval( refmat_ ), R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( repeat( eval( omat_ )  , R0, R1 ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( repeat( eval( refmat_ ), R0, R1 ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Submatrix-wise repeat with Schur product assignment with evaluated matrix (compile time)
      {
         test_  = "Submatrix-wise repeat with Schur product assignment with evaluated matrix (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( repeat<R0,R1>( eval( mat_ )    ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( repeat<R0,R1>( eval( refmat_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<mat_.rows()*R0; row+=m ) {
               m = blaze::rand<size_t>( 1UL, mat_.rows()*R0 - row );
               for( size_t column=0UL, n=0UL; column<mat_.columns()*R1; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, mat_.columns()*R1 - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( repeat<R0,R1>( eval( omat_ )   ), row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( repeat<R0,R1>( eval( refmat_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds submatrix construction (invalid number of rows)
      {
         test_  = "Out-of-bounds submatrix construction (invalid number of rows)";
         error_ = "Setup of out-of-bounds submatrix succeeded";

         try {
            auto sm = submatrix( repeat( mat_, R0, R1 ), 1UL, 0UL, mat_.rows()*R0, mat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( repeat<R0,R1>( mat_ ), 1UL, 0UL, mat_.rows()*R0, mat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( repeat( omat_, R0, R1 ), 1UL, 0UL, omat_.rows()*R0, omat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( repeat<R0,R1>( omat_ ), 1UL, 0UL, omat_.rows()*R0, omat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
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
            auto sm = submatrix( repeat( mat_, R0, R1 ), 0UL, 1UL, mat_.rows()*R0, mat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( repeat<R0,R1>( mat_ ), 0UL, 1UL, mat_.rows()*R0, mat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( repeat( omat_, R0, R1 ), 0UL, 1UL, omat_.rows()*R0, omat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( repeat<R0,R1>( omat_ ), 0UL, 1UL, omat_.rows()*R0, omat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
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
            auto sm = submatrix( repeat( mat_, R0, R1 ), mat_.rows()*R0, 0UL, 1UL, mat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( repeat<R0,R1>( mat_ ), mat_.rows()*R0, 0UL, 1UL, mat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( repeat( omat_, R0, R1 ), omat_.rows()*R0, 0UL, 1UL, omat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( repeat<R0,R1>( omat_ ), omat_.rows()*R0, 0UL, 1UL, omat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
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
            auto sm = submatrix( repeat( mat_, R0, R1 ), 0UL, mat_.columns()*R1, mat_.rows()*R0, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( repeat<R0,R1>( mat_ ), 0UL, mat_.columns()*R1, mat_.rows()*R0, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( repeat( omat_, R0, R1 ), 0UL, omat_.columns()*R1, omat_.rows()*R0, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( repeat<R0,R1>( omat_ ), 0UL, omat_.columns()*R1, omat_.rows()*R0, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
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
/*!\brief Testing the row-wise sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the row-wise matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the repeat operation or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testRowOperation()
{
#if BLAZETEST_MATHTEST_TEST_ROW_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ROW_OPERATION > 1 )
   {
      using blaze::repeat;


      if( mat_.rows() == 0UL )
         return;


      //=====================================================================================
      // Row-wise repeat operation
      //=====================================================================================

      // Row-wise repeat operation with the given matrix (runtime)
      {
         test_  = "Row-wise repeat operation with the given matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) = row( repeat( mat_   , R0, R1 ), i );
               row( odres_ , i ) = row( repeat( mat_   , R0, R1 ), i );
               row( sres_  , i ) = row( repeat( mat_   , R0, R1 ), i );
               row( osres_ , i ) = row( repeat( mat_   , R0, R1 ), i );
               row( refres_, i ) = row( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) = row( repeat( omat_  , R0, R1 ), i );
               row( odres_ , i ) = row( repeat( omat_  , R0, R1 ), i );
               row( sres_  , i ) = row( repeat( omat_  , R0, R1 ), i );
               row( osres_ , i ) = row( repeat( omat_  , R0, R1 ), i );
               row( refres_, i ) = row( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Row-wise repeat operation with the given matrix (compile time)
      {
         test_  = "Row-wise repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) = row( repeat<R0,R1>( mat_    ), i );
               row( odres_ , i ) = row( repeat<R0,R1>( mat_    ), i );
               row( sres_  , i ) = row( repeat<R0,R1>( mat_    ), i );
               row( osres_ , i ) = row( repeat<R0,R1>( mat_    ), i );
               row( refres_, i ) = row( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) = row( repeat<R0,R1>( omat_   ), i );
               row( odres_ , i ) = row( repeat<R0,R1>( omat_   ), i );
               row( sres_  , i ) = row( repeat<R0,R1>( omat_   ), i );
               row( osres_ , i ) = row( repeat<R0,R1>( omat_   ), i );
               row( refres_, i ) = row( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Row-wise repeat operation with evaluated matrix (runtime)
      {
         test_  = "Row-wise repeat operation with evaluated matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) = row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( odres_ , i ) = row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( sres_  , i ) = row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( osres_ , i ) = row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( refres_, i ) = row( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) = row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( odres_ , i ) = row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( sres_  , i ) = row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( osres_ , i ) = row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( refres_, i ) = row( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Row-wise repeat operation with evaluated matrix (compile time)
      {
         test_  = "Row-wise repeat operation with evaluated matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) = row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( odres_ , i ) = row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( sres_  , i ) = row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( osres_ , i ) = row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( refres_, i ) = row( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) = row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( odres_ , i ) = row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( sres_  , i ) = row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( osres_ , i ) = row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( refres_, i ) = row( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Row-wise repeat with addition assignment
      //=====================================================================================

      // Row-wise repeat with addition assignment with the given matrix (runtime)
      {
         test_  = "Row-wise repeat with addition assignment with the given matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) += row( repeat( mat_   , R0, R1 ), i );
               row( odres_ , i ) += row( repeat( mat_   , R0, R1 ), i );
               row( sres_  , i ) += row( repeat( mat_   , R0, R1 ), i );
               row( osres_ , i ) += row( repeat( mat_   , R0, R1 ), i );
               row( refres_, i ) += row( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) += row( repeat( omat_  , R0, R1 ), i );
               row( odres_ , i ) += row( repeat( omat_  , R0, R1 ), i );
               row( sres_  , i ) += row( repeat( omat_  , R0, R1 ), i );
               row( osres_ , i ) += row( repeat( omat_  , R0, R1 ), i );
               row( refres_, i ) += row( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Row-wise repeat with addition assignment with the given matrix (compile time)
      {
         test_  = "Row-wise repeat with addition assignment with the given matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) += row( repeat<R0,R1>( mat_    ), i );
               row( odres_ , i ) += row( repeat<R0,R1>( mat_    ), i );
               row( sres_  , i ) += row( repeat<R0,R1>( mat_    ), i );
               row( osres_ , i ) += row( repeat<R0,R1>( mat_    ), i );
               row( refres_, i ) += row( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) += row( repeat<R0,R1>( omat_   ), i );
               row( odres_ , i ) += row( repeat<R0,R1>( omat_   ), i );
               row( sres_  , i ) += row( repeat<R0,R1>( omat_   ), i );
               row( osres_ , i ) += row( repeat<R0,R1>( omat_   ), i );
               row( refres_, i ) += row( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Row-wise repeat with addition assignment with evaluated matrix (runtime)
      {
         test_  = "Row-wise repeat with addition assignment with evaluated matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) += row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( odres_ , i ) += row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( sres_  , i ) += row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( osres_ , i ) += row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( refres_, i ) += row( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) += row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( odres_ , i ) += row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( sres_  , i ) += row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( osres_ , i ) += row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( refres_, i ) += row( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Row-wise repeat with addition assignment with evaluated matrix (compile time)
      {
         test_  = "Row-wise repeat with addition assignment with evaluated matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) += row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( odres_ , i ) += row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( sres_  , i ) += row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( osres_ , i ) += row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( refres_, i ) += row( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) += row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( odres_ , i ) += row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( sres_  , i ) += row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( osres_ , i ) += row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( refres_, i ) += row( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Row-wise repeat with subtraction assignment
      //=====================================================================================

      // Row-wise repeat with subtraction assignment with the given matrix (runtime)
      {
         test_  = "Row-wise repeat with subtraction assignment with the given matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) -= row( repeat( mat_   , R0, R1 ), i );
               row( odres_ , i ) -= row( repeat( mat_   , R0, R1 ), i );
               row( sres_  , i ) -= row( repeat( mat_   , R0, R1 ), i );
               row( osres_ , i ) -= row( repeat( mat_   , R0, R1 ), i );
               row( refres_, i ) -= row( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) -= row( repeat( omat_  , R0, R1 ), i );
               row( odres_ , i ) -= row( repeat( omat_  , R0, R1 ), i );
               row( sres_  , i ) -= row( repeat( omat_  , R0, R1 ), i );
               row( osres_ , i ) -= row( repeat( omat_  , R0, R1 ), i );
               row( refres_, i ) -= row( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Row-wise repeat with subtraction assignment with the given matrix (compile time)
      {
         test_  = "Row-wise repeat with subtraction assignment with the given matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) -= row( repeat<R0,R1>( mat_    ), i );
               row( odres_ , i ) -= row( repeat<R0,R1>( mat_    ), i );
               row( sres_  , i ) -= row( repeat<R0,R1>( mat_    ), i );
               row( osres_ , i ) -= row( repeat<R0,R1>( mat_    ), i );
               row( refres_, i ) -= row( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) -= row( repeat<R0,R1>( omat_   ), i );
               row( odres_ , i ) -= row( repeat<R0,R1>( omat_   ), i );
               row( sres_  , i ) -= row( repeat<R0,R1>( omat_   ), i );
               row( osres_ , i ) -= row( repeat<R0,R1>( omat_   ), i );
               row( refres_, i ) -= row( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Row-wise repeat with subtraction assignment with evaluated matrix (runtime)
      {
         test_  = "Row-wise repeat with subtraction assignment with evaluated matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) -= row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( odres_ , i ) -= row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( sres_  , i ) -= row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( osres_ , i ) -= row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( refres_, i ) -= row( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) -= row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( odres_ , i ) -= row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( sres_  , i ) -= row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( osres_ , i ) -= row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( refres_, i ) -= row( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Row-wise repeat with subtraction assignment with evaluated matrix (compile time)
      {
         test_  = "Row-wise repeat with subtraction assignment with evaluated matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) -= row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( odres_ , i ) -= row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( sres_  , i ) -= row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( osres_ , i ) -= row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( refres_, i ) -= row( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) -= row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( odres_ , i ) -= row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( sres_  , i ) -= row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( osres_ , i ) -= row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( refres_, i ) -= row( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Row-wise repeat with multiplication assignment
      //=====================================================================================

      // Row-wise repeat with multiplication assignment with the given matrix (runtime)
      {
         test_  = "Row-wise repeat with multiplication assignment with the given matrix (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) *= row( repeat( mat_   , R0, R1 ), i );
               row( odres_ , i ) *= row( repeat( mat_   , R0, R1 ), i );
               row( sres_  , i ) *= row( repeat( mat_   , R0, R1 ), i );
               row( osres_ , i ) *= row( repeat( mat_   , R0, R1 ), i );
               row( refres_, i ) *= row( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) *= row( repeat( omat_  , R0, R1 ), i );
               row( odres_ , i ) *= row( repeat( omat_  , R0, R1 ), i );
               row( sres_  , i ) *= row( repeat( omat_  , R0, R1 ), i );
               row( osres_ , i ) *= row( repeat( omat_  , R0, R1 ), i );
               row( refres_, i ) *= row( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Row-wise repeat with multiplication assignment with the given matrix (compile time)
      {
         test_  = "Row-wise repeat with multiplication assignment with the given matrix (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) *= row( repeat<R0,R1>( mat_    ), i );
               row( odres_ , i ) *= row( repeat<R0,R1>( mat_    ), i );
               row( sres_  , i ) *= row( repeat<R0,R1>( mat_    ), i );
               row( osres_ , i ) *= row( repeat<R0,R1>( mat_    ), i );
               row( refres_, i ) *= row( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) *= row( repeat<R0,R1>( omat_   ), i );
               row( odres_ , i ) *= row( repeat<R0,R1>( omat_   ), i );
               row( sres_  , i ) *= row( repeat<R0,R1>( omat_   ), i );
               row( osres_ , i ) *= row( repeat<R0,R1>( omat_   ), i );
               row( refres_, i ) *= row( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Row-wise repeat with multiplication assignment with evaluated matrix (runtime)
      {
         test_  = "Row-wise repeat with multiplication assignment with evaluated matrix (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) *= row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( odres_ , i ) *= row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( sres_  , i ) *= row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( osres_ , i ) *= row( repeat( eval( mat_ )   , R0, R1 ), i );
               row( refres_, i ) *= row( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) *= row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( odres_ , i ) *= row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( sres_  , i ) *= row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( osres_ , i ) *= row( repeat( eval( omat_ )  , R0, R1 ), i );
               row( refres_, i ) *= row( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Row-wise repeat with multiplication assignment with evaluated matrix (compile time)
      {
         test_  = "Row-wise repeat with multiplication assignment with evaluated matrix (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) *= row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( odres_ , i ) *= row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( sres_  , i ) *= row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( osres_ , i ) *= row( repeat<R0,R1>( eval( mat_ )    ), i );
               row( refres_, i ) *= row( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t i=0UL; i<mat_.rows(); ++i ) {
               row( dres_  , i ) *= row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( odres_ , i ) *= row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( sres_  , i ) *= row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( osres_ , i ) *= row( repeat<R0,R1>( eval( omat_ )   ), i );
               row( refres_, i ) *= row( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid row index)
      {
         test_  = "Out-of-bounds row construction (invalid row index)";
         error_ = "Setup of out-of-bounds row succeeded";

         try {
            auto r = row( repeat( mat_, R0, R1 ), mat_.rows()*R0 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = row( repeat<R0,R1>( mat_ ), mat_.rows()*R0 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = row( repeat( omat_, R0, R1 ), omat_.rows()*R0 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = row( repeat<R0,R1>( omat_ ), omat_.rows()*R0 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
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
/*!\brief Testing the rows-wise sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the rows-wise matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the repeat operation or the subsequent assignment is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testRowsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_ROWS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ROWS_OPERATION > 1 )
   {
      using blaze::repeat;


      if( mat_.rows() == 0UL )
         return;


      std::vector<size_t> indices( mat_.rows() );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );


      //=====================================================================================
      // Rows-wise repeat operation
      //=====================================================================================

      // Rows-wise repeat operation with the given matrix (runtime)
      {
         test_  = "Rows-wise repeat operation with the given matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( repeat( refmat_, R0, R1 ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) = rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( repeat( refmat_, R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Rows-wise repeat operation with the given matrix (compile time)
      {
         test_  = "Rows-wise repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( repeat<R0,R1>( refmat_ ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) = rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( repeat<R0,R1>( refmat_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Rows-wise repeat operation with evaluated matrix (runtime)
      {
         test_  = "Rows-wise repeat operation with evaluated matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) = rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Rows-wise repeat operation with evaluated matrix (compile time)
      {
         test_  = "Rows-wise repeat operation with evaluated matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) = rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Rows-wise repeat with addition assignment
      //=====================================================================================

      // Rows-wise repeat with addition assignment with the given matrix (runtime)
      {
         test_  = "Rows-wise repeat with addition assignment with the given matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( repeat( refmat_, R0, R1 ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) += rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( repeat( refmat_, R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Rows-wise repeat with addition assignment with the given matrix (compile time)
      {
         test_  = "Rows-wise repeat with addition assignment with the given matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( repeat<R0,R1>( refmat_ ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) += rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( repeat<R0,R1>( refmat_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Rows-wise repeat with addition assignment with evaluated matrix (runtime)
      {
         test_  = "Rows-wise repeat with addition assignment with evaluated matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) += rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Rows-wise repeat with addition assignment with evaluated matrix (compile time)
      {
         test_  = "Rows-wise repeat with addition assignment with evaluated matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) += rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Rows-wise repeat with subtraction assignment
      //=====================================================================================

      // Rows-wise repeat with subtraction assignment with the given matrix (runtime)
      {
         test_  = "Rows-wise repeat with subtraction assignment with the given matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( repeat( refmat_, R0, R1 ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) -= rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( repeat( refmat_, R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Rows-wise repeat with subtraction assignment with the given matrix (compile time)
      {
         test_  = "Rows-wise repeat with subtraction assignment with the given matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( repeat<R0,R1>( refmat_ ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) -= rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( repeat<R0,R1>( refmat_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Rows-wise repeat with subtraction assignment with evaluated matrix (runtime)
      {
         test_  = "Rows-wise repeat with subtraction assignment with evaluated matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) -= rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Rows-wise repeat with subtraction assignment with evaluated matrix (compile time)
      {
         test_  = "Rows-wise repeat with subtraction assignment with evaluated matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) -= rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Rows-wise repeat with Schur product assignment
      //=====================================================================================

      // Rows-wise repeat with Schur product assignment with the given matrix (runtime)
      {
         test_  = "Rows-wise repeat with Schur product assignment with the given matrix (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( repeat( mat_   , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( repeat( refmat_, R0, R1 ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) %= rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( repeat( omat_  , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( repeat( refmat_, R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Rows-wise repeat with Schur product assignment with the given matrix (compile time)
      {
         test_  = "Rows-wise repeat with Schur product assignment with the given matrix (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( repeat<R0,R1>( mat_    ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( repeat<R0,R1>( refmat_ ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) %= rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( repeat<R0,R1>( omat_   ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( repeat<R0,R1>( refmat_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Rows-wise repeat with Schur product assignment with evaluated matrix (runtime)
      {
         test_  = "Rows-wise repeat with Schur product assignment with evaluated matrix (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) %= rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Rows-wise repeat with Schur product assignment with evaluated matrix (compile time)
      {
         test_  = "Rows-wise repeat with Schur product assignment with evaluated matrix (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
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
               rows( dres_  , &indices[index], n ) %= rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid row index; initializer_list)
      {
         test_  = "Out-of-bounds row selection construction (invalid row index; initializer_list)";
         error_ = "Setup of out-of-bounds row selection succeeded";

         try {
            auto r = rows( repeat( mat_, R0, R1 ), { mat_.rows()*R0 } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( repeat<R0,R1>( mat_ ), { mat_.rows()*R0 } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( repeat( omat_, R0, R1 ), { omat_.rows()*R0 } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( repeat<R0,R1>( omat_ ), { omat_.rows()*R0 } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
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
            auto r = rows( repeat( mat_, R0, R1 ), [index=mat_.rows()*R0]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( repeat<R0,R1>( mat_ ), [index=mat_.rows()*R0]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( repeat( omat_, R0, R1 ), [index=omat_.rows()*R0]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( repeat<R0,R1>( omat_ ), [index=omat_.rows()*R0]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
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
/*!\brief Skipping the rows-wise sparse matrix repeat operation.
//
// \return void
//
// This function is called in case the rows-wise matrix repeat operation is not available for
// the given matrix type \a MT.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testRowsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the column-wise sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the column-wise matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the addition or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testColumnOperation()
{
#if BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION
   if( BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION > 1 )
   {
      using blaze::repeat;


      if( mat_.columns() == 0UL )
         return;


      //=====================================================================================
      // Column-wise repeat operation
      //=====================================================================================

      // Column-wise repeat operation with the given matrix (runtime)
      {
         test_  = "Column-wise repeat operation with the given matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) = column( repeat( mat_   , R0, R1 ), j );
               column( odres_ , j ) = column( repeat( mat_   , R0, R1 ), j );
               column( sres_  , j ) = column( repeat( mat_   , R0, R1 ), j );
               column( osres_ , j ) = column( repeat( mat_   , R0, R1 ), j );
               column( refres_, j ) = column( repeat( refmat_, R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) = column( repeat( omat_  , R0, R1 ), j );
               column( odres_ , j ) = column( repeat( omat_  , R0, R1 ), j );
               column( sres_  , j ) = column( repeat( omat_  , R0, R1 ), j );
               column( osres_ , j ) = column( repeat( omat_  , R0, R1 ), j );
               column( refres_, j ) = column( repeat( refmat_, R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Column-wise repeat operation with the given matrix (compile time)
      {
         test_  = "Column-wise repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) = column( repeat<R0,R1>( mat_    ), j );
               column( odres_ , j ) = column( repeat<R0,R1>( mat_    ), j );
               column( sres_  , j ) = column( repeat<R0,R1>( mat_    ), j );
               column( osres_ , j ) = column( repeat<R0,R1>( mat_    ), j );
               column( refres_, j ) = column( repeat<R0,R1>( refmat_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) = column( repeat<R0,R1>( omat_   ), j );
               column( odres_ , j ) = column( repeat<R0,R1>( omat_   ), j );
               column( sres_  , j ) = column( repeat<R0,R1>( omat_   ), j );
               column( osres_ , j ) = column( repeat<R0,R1>( omat_   ), j );
               column( refres_, j ) = column( repeat<R0,R1>( refmat_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Column-wise repeat operation with evaluated matrix (runtime)
      {
         test_  = "Column-wise repeat operation with evaluated matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) = column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( odres_ , j ) = column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( sres_  , j ) = column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( osres_ , j ) = column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( refres_, j ) = column( repeat( eval( refmat_ ), R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) = column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( odres_ , j ) = column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( sres_  , j ) = column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( osres_ , j ) = column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( refres_, j ) = column( repeat( eval( refmat_ ), R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Column-wise repeat operation with evaluated matrix (compile time)
      {
         test_  = "Column-wise repeat operation with evaluated matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) = column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( odres_ , j ) = column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( sres_  , j ) = column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( osres_ , j ) = column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( refres_, j ) = column( repeat<R0,R1>( eval( refmat_ ) ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) = column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( odres_ , j ) = column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( sres_  , j ) = column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( osres_ , j ) = column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( refres_, j ) = column( repeat<R0,R1>( eval( refmat_ ) ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Column-wise repeat with addition assignment
      //=====================================================================================

      // Column-wise repeat with addition assignment with the given matrix (runtime)
      {
         test_  = "Column-wise repeat with addition assignment with the given matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) += column( repeat( mat_   , R0, R1 ), j );
               column( odres_ , j ) += column( repeat( mat_   , R0, R1 ), j );
               column( sres_  , j ) += column( repeat( mat_   , R0, R1 ), j );
               column( osres_ , j ) += column( repeat( mat_   , R0, R1 ), j );
               column( refres_, j ) += column( repeat( refmat_, R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) += column( repeat( omat_  , R0, R1 ), j );
               column( odres_ , j ) += column( repeat( omat_  , R0, R1 ), j );
               column( sres_  , j ) += column( repeat( omat_  , R0, R1 ), j );
               column( osres_ , j ) += column( repeat( omat_  , R0, R1 ), j );
               column( refres_, j ) += column( repeat( refmat_, R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Column-wise repeat with addition assignment with the given matrix (compile time)
      {
         test_  = "Column-wise repeat with addition assignment with the given matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) += column( repeat<R0,R1>( mat_    ), j );
               column( odres_ , j ) += column( repeat<R0,R1>( mat_    ), j );
               column( sres_  , j ) += column( repeat<R0,R1>( mat_    ), j );
               column( osres_ , j ) += column( repeat<R0,R1>( mat_    ), j );
               column( refres_, j ) += column( repeat<R0,R1>( refmat_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) += column( repeat<R0,R1>( omat_   ), j );
               column( odres_ , j ) += column( repeat<R0,R1>( omat_   ), j );
               column( sres_  , j ) += column( repeat<R0,R1>( omat_   ), j );
               column( osres_ , j ) += column( repeat<R0,R1>( omat_   ), j );
               column( refres_, j ) += column( repeat<R0,R1>( refmat_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Column-wise repeat with addition assignment with evaluated matrix (runtime)
      {
         test_  = "Column-wise repeat with addition assignment with evaluated matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) += column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( odres_ , j ) += column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( sres_  , j ) += column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( osres_ , j ) += column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( refres_, j ) += column( repeat( eval( refmat_ ), R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) += column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( odres_ , j ) += column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( sres_  , j ) += column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( osres_ , j ) += column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( refres_, j ) += column( repeat( eval( refmat_ ), R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Column-wise repeat with addition assignment with evaluated matrix (compile time)
      {
         test_  = "Column-wise repeat with addition assignment with evaluated matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) += column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( odres_ , j ) += column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( sres_  , j ) += column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( osres_ , j ) += column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( refres_, j ) += column( repeat<R0,R1>( eval( refmat_ ) ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) += column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( odres_ , j ) += column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( sres_  , j ) += column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( osres_ , j ) += column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( refres_, j ) += column( repeat<R0,R1>( eval( refmat_ ) ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Column-wise repeat with subtraction assignment
      //=====================================================================================

      // Column-wise repeat with subtraction assignment with the given matrix (runtime)
      {
         test_  = "Column-wise repeat with subtraction assignment with the given matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat( mat_   , R0, R1 ), j );
               column( odres_ , j ) -= column( repeat( mat_   , R0, R1 ), j );
               column( sres_  , j ) -= column( repeat( mat_   , R0, R1 ), j );
               column( osres_ , j ) -= column( repeat( mat_   , R0, R1 ), j );
               column( refres_, j ) -= column( repeat( refmat_, R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat( omat_  , R0, R1 ), j );
               column( odres_ , j ) -= column( repeat( omat_  , R0, R1 ), j );
               column( sres_  , j ) -= column( repeat( omat_  , R0, R1 ), j );
               column( osres_ , j ) -= column( repeat( omat_  , R0, R1 ), j );
               column( refres_, j ) -= column( repeat( refmat_, R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Column-wise repeat with subtraction assignment with the given matrix (compile time)
      {
         test_  = "Column-wise repeat with subtraction assignment with the given matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat<R0,R1>( mat_    ), j );
               column( odres_ , j ) -= column( repeat<R0,R1>( mat_    ), j );
               column( sres_  , j ) -= column( repeat<R0,R1>( mat_    ), j );
               column( osres_ , j ) -= column( repeat<R0,R1>( mat_    ), j );
               column( refres_, j ) -= column( repeat<R0,R1>( refmat_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat<R0,R1>( omat_   ), j );
               column( odres_ , j ) -= column( repeat<R0,R1>( omat_   ), j );
               column( sres_  , j ) -= column( repeat<R0,R1>( omat_   ), j );
               column( osres_ , j ) -= column( repeat<R0,R1>( omat_   ), j );
               column( refres_, j ) -= column( repeat<R0,R1>( refmat_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Column-wise repeat with subtraction assignment with evaluated matrix (runtime)
      {
         test_  = "Column-wise repeat with subtraction assignment with evaluated matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( odres_ , j ) -= column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( sres_  , j ) -= column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( osres_ , j ) -= column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( refres_, j ) -= column( repeat( eval( refmat_ ), R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( odres_ , j ) -= column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( sres_  , j ) -= column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( osres_ , j ) -= column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( refres_, j ) -= column( repeat( eval( refmat_ ), R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Column-wise repeat with subtraction assignment with evaluated matrix (compile time)
      {
         test_  = "Column-wise repeat with subtraction assignment with evaluated matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( odres_ , j ) -= column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( sres_  , j ) -= column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( osres_ , j ) -= column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( refres_, j ) -= column( repeat<R0,R1>( eval( refmat_ ) ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( odres_ , j ) -= column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( sres_  , j ) -= column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( osres_ , j ) -= column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( refres_, j ) -= column( repeat<R0,R1>( eval( refmat_ ) ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Column-wise repeat with multiplication assignment
      //=====================================================================================

      // Column-wise repeat with multiplication assignment with the given matrix (runtime)
      {
         test_  = "Column-wise repeat with multiplication assignment with the given matrix (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat( mat_   , R0, R1 ), j );
               column( odres_ , j ) -= column( repeat( mat_   , R0, R1 ), j );
               column( sres_  , j ) -= column( repeat( mat_   , R0, R1 ), j );
               column( osres_ , j ) -= column( repeat( mat_   , R0, R1 ), j );
               column( refres_, j ) -= column( repeat( refmat_, R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat( omat_  , R0, R1 ), j );
               column( odres_ , j ) -= column( repeat( omat_  , R0, R1 ), j );
               column( sres_  , j ) -= column( repeat( omat_  , R0, R1 ), j );
               column( osres_ , j ) -= column( repeat( omat_  , R0, R1 ), j );
               column( refres_, j ) -= column( repeat( refmat_, R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Column-wise repeat with multiplication assignment with the given matrix (compile time)
      {
         test_  = "Column-wise repeat with multiplication assignment with the given matrix (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat<R0,R1>( mat_    ), j );
               column( odres_ , j ) -= column( repeat<R0,R1>( mat_    ), j );
               column( sres_  , j ) -= column( repeat<R0,R1>( mat_    ), j );
               column( osres_ , j ) -= column( repeat<R0,R1>( mat_    ), j );
               column( refres_, j ) -= column( repeat<R0,R1>( refmat_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat<R0,R1>( omat_   ), j );
               column( odres_ , j ) -= column( repeat<R0,R1>( omat_   ), j );
               column( sres_  , j ) -= column( repeat<R0,R1>( omat_   ), j );
               column( osres_ , j ) -= column( repeat<R0,R1>( omat_   ), j );
               column( refres_, j ) -= column( repeat<R0,R1>( refmat_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Column-wise repeat with multiplication assignment with evaluated matrix (runtime)
      {
         test_  = "Column-wise repeat with multiplication assignment with evaluated matrix (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( odres_ , j ) -= column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( sres_  , j ) -= column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( osres_ , j ) -= column( repeat( eval( mat_ )   , R0, R1 ), j );
               column( refres_, j ) -= column( repeat( eval( refmat_ ), R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( odres_ , j ) -= column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( sres_  , j ) -= column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( osres_ , j ) -= column( repeat( eval( omat_ )  , R0, R1 ), j );
               column( refres_, j ) -= column( repeat( eval( refmat_ ), R0, R1 ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Column-wise repeat with multiplication assignment with evaluated matrix (compile time)
      {
         test_  = "Column-wise repeat with multiplication assignment with evaluated matrix (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( odres_ , j ) -= column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( sres_  , j ) -= column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( osres_ , j ) -= column( repeat<R0,R1>( eval( mat_ )    ), j );
               column( refres_, j ) -= column( repeat<R0,R1>( eval( refmat_ ) ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( size_t j=0UL; j<mat_.columns(); ++j ) {
               column( dres_  , j ) -= column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( odres_ , j ) -= column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( sres_  , j ) -= column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( osres_ , j ) -= column( repeat<R0,R1>( eval( omat_ )   ), j );
               column( refres_, j ) -= column( repeat<R0,R1>( eval( refmat_ ) ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid column index)
      {
         test_  = "Out-of-bounds column construction (invalid column index)";
         error_ = "Setup of out-of-bounds column succeeded";

         try {
            auto c = column( repeat( mat_, R0, R1 ), mat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = column( repeat<R0,R1>( mat_ ), mat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = column( repeat( omat_, R0, R1 ), omat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = column( repeat<R0,R1>( omat_ ), omat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
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
/*!\brief Testing the columns-wise sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the columns-wise matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the addition or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testColumnsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION > 1 )
   {
      using blaze::repeat;


      if( mat_.columns() == 0UL )
         return;


      std::vector<size_t> indices( mat_.columns() );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );


      //=====================================================================================
      // Columns-wise repeat operation
      //=====================================================================================

      // Columns-wise repeat operation with the given matrix (runtime)
      {
         test_  = "Columns-wise repeat operation with the given matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( repeat( refmat_, R0, R1 ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) = columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( repeat( refmat_, R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Columns-wise repeat operation with the given matrix (compile time)
      {
         test_  = "Columns-wise repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( repeat<R0,R1>( refmat_ ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) = columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( repeat<R0,R1>( refmat_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Columns-wise repeat operation with evaluated matrix (runtime)
      {
         test_  = "Columns-wise repeat operation with evaluated matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) = columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Columns-wise repeat operation with evaluated matrix (compile time)
      {
         test_  = "Columns-wise repeat operation with evaluated matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) = columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Columns-wise repeat with addition assignment
      //=====================================================================================

      // Columns-wise repeat with addition assignment with the given matrix (runtime)
      {
         test_  = "Columns-wise repeat with addition assignment with the given matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( repeat( refmat_, R0, R1 ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) += columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( repeat( refmat_, R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Columns-wise repeat with addition assignment with the given matrix (compile time)
      {
         test_  = "Columns-wise repeat with addition assignment with the given matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( repeat<R0,R1>( refmat_ ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) += columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( repeat<R0,R1>( refmat_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Columns-wise repeat with addition assignment with evaluated matrix (runtime)
      {
         test_  = "Columns-wise repeat with addition assignment with evaluated matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) += columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Columns-wise repeat with addition assignment with evaluated matrix (compile time)
      {
         test_  = "Columns-wise repeat with addition assignment with evaluated matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) += columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Columns-wise repeat with subtraction assignment
      //=====================================================================================

      // Columns-wise repeat with subtraction assignment with the given matrix (runtime)
      {
         test_  = "Columns-wise repeat with subtraction assignment with the given matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( repeat( refmat_, R0, R1 ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) -= columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( repeat( refmat_, R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Columns-wise repeat with subtraction assignment with the given matrix (compile time)
      {
         test_  = "Columns-wise repeat with subtraction assignment with the given matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( repeat<R0,R1>( refmat_ ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) -= columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( repeat<R0,R1>( refmat_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Columns-wise repeat with subtraction assignment with evaluated matrix (runtime)
      {
         test_  = "Columns-wise repeat with subtraction assignment with evaluated matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) -= columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Columns-wise repeat with subtraction assignment with evaluated matrix (compile time)
      {
         test_  = "Columns-wise repeat with subtraction assignment with evaluated matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) -= columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Columns-wise repeat with Schur product assignment
      //=====================================================================================

      // Columns-wise repeat with Schur product assignment with the given matrix (runtime)
      {
         test_  = "Columns-wise repeat with Schur product assignment with the given matrix (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( repeat( mat_   , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( repeat( refmat_, R0, R1 ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) %= columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( repeat( omat_  , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( repeat( refmat_, R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Columns-wise repeat with Schur product assignment with the given matrix (compile time)
      {
         test_  = "Columns-wise repeat with Schur product assignment with the given matrix (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( repeat<R0,R1>( mat_    ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( repeat<R0,R1>( refmat_ ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) %= columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( repeat<R0,R1>( omat_   ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( repeat<R0,R1>( refmat_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Columns-wise repeat with Schur product assignment with evaluated matrix (runtime)
      {
         test_  = "Columns-wise repeat with Schur product assignment with evaluated matrix (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( repeat( eval( mat_ )   , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) %= columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( repeat( eval( omat_ )  , R0, R1 ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( repeat( eval( refmat_ ), R0, R1 ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Columns-wise repeat with Schur product assignment with evaluated matrix (compile time)
      {
         test_  = "Columns-wise repeat with Schur product assignment with evaluated matrix (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( repeat<R0,R1>( eval( mat_ )    ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
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
               columns( dres_  , &indices[index], n ) %= columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( repeat<R0,R1>( eval( omat_ )   ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( repeat<R0,R1>( eval( refmat_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid column index; initializer_list)
      {
         test_  = "Out-of-bounds column selection construction (invalid column index; initializer_list)";
         error_ = "Setup of out-of-bounds column selection succeeded";

         try {
            auto c = columns( repeat( mat_, R0, R1 ), { mat_.columns()*R1 } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( repeat<R0,R1>( mat_ ), { mat_.columns()*R1 } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( repeat( omat_, R0, R1 ), { omat_.columns()*R1 } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( repeat<R0,R1>( omat_ ), { omat_.columns()*R1 } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
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
            auto c = columns( repeat( mat_, R0, R1 ), [index=mat_.columns()*R1]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( repeat<R0,R1>( mat_ ), [index=mat_.columns()*R1]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( repeat( omat_, R0, R1 ), [index=omat_.columns()*R1]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( repeat<R0,R1>( omat_ ), [index=omat_.columns()*R1]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
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
/*!\brief Skipping the columns-wise sparse matrix repeat operation.
//
// \return void
//
// This function is called in case the columns-wise matrix repeat operation is not available for
// the given matrix type \a MT.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testColumnsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the band-wise sparse matrix repeat operation.
//
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the band-wise matrix repeat operation with plain assignment, addition
// assignment, subtraction assignment, and multiplication assignment. In case any error resulting
// from the addition or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::testBandOperation()
{
#if BLAZETEST_MATHTEST_TEST_BAND_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BAND_OPERATION > 1 )
   {
      using blaze::repeat;


      if( mat_.rows() == 0UL || mat_.columns() == 0UL )
         return;


      const ptrdiff_t ibegin( 1UL - mat_.rows() );
      const ptrdiff_t iend  ( mat_.columns() );


      //=====================================================================================
      // Band-wise repeat operation
      //=====================================================================================

      // Band-wise repeat operation with the given matrix (runtime)
      {
         test_  = "Band-wise repeat operation with the given matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( repeat( mat_   , R0, R1 ), i );
               band( odres_ , i ) = band( repeat( mat_   , R0, R1 ), i );
               band( sres_  , i ) = band( repeat( mat_   , R0, R1 ), i );
               band( osres_ , i ) = band( repeat( mat_   , R0, R1 ), i );
               band( refres_, i ) = band( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( repeat( omat_  , R0, R1 ), i );
               band( odres_ , i ) = band( repeat( omat_  , R0, R1 ), i );
               band( sres_  , i ) = band( repeat( omat_  , R0, R1 ), i );
               band( osres_ , i ) = band( repeat( omat_  , R0, R1 ), i );
               band( refres_, i ) = band( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Band-wise repeat operation with the given matrix (compile time)
      {
         test_  = "Band-wise repeat operation with the given matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( repeat<R0,R1>( mat_    ), i );
               band( odres_ , i ) = band( repeat<R0,R1>( mat_    ), i );
               band( sres_  , i ) = band( repeat<R0,R1>( mat_    ), i );
               band( osres_ , i ) = band( repeat<R0,R1>( mat_    ), i );
               band( refres_, i ) = band( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( repeat<R0,R1>( omat_   ), i );
               band( odres_ , i ) = band( repeat<R0,R1>( omat_   ), i );
               band( sres_  , i ) = band( repeat<R0,R1>( omat_   ), i );
               band( osres_ , i ) = band( repeat<R0,R1>( omat_   ), i );
               band( refres_, i ) = band( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Band-wise repeat operation with evaluated matrix (runtime)
      {
         test_  = "Band-wise repeat operation with evaluated matrix (runtime)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( odres_ , i ) = band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( sres_  , i ) = band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( osres_ , i ) = band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( refres_, i ) = band( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( odres_ , i ) = band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( sres_  , i ) = band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( osres_ , i ) = band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( refres_, i ) = band( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Band-wise repeat operation with evaluated matrix (compile time)
      {
         test_  = "Band-wise repeat operation with evaluated matrix (compile time)";
         error_ = "Failed repeat operation";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( odres_ , i ) = band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( sres_  , i ) = band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( osres_ , i ) = band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( refres_, i ) = band( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) = band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( odres_ , i ) = band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( sres_  , i ) = band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( osres_ , i ) = band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( refres_, i ) = band( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Band-wise repeat with addition assignment
      //=====================================================================================

      // Band-wise repeat with addition assignment with the given matrix (runtime)
      {
         test_  = "Band-wise repeat with addition assignment with the given matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( repeat( mat_   , R0, R1 ), i );
               band( odres_ , i ) += band( repeat( mat_   , R0, R1 ), i );
               band( sres_  , i ) += band( repeat( mat_   , R0, R1 ), i );
               band( osres_ , i ) += band( repeat( mat_   , R0, R1 ), i );
               band( refres_, i ) += band( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( repeat( omat_  , R0, R1 ), i );
               band( odres_ , i ) += band( repeat( omat_  , R0, R1 ), i );
               band( sres_  , i ) += band( repeat( omat_  , R0, R1 ), i );
               band( osres_ , i ) += band( repeat( omat_  , R0, R1 ), i );
               band( refres_, i ) += band( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Band-wise repeat with addition assignment with the given matrix (compile time)
      {
         test_  = "Band-wise repeat with addition assignment with the given matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( repeat<R0,R1>( mat_    ), i );
               band( odres_ , i ) += band( repeat<R0,R1>( mat_    ), i );
               band( sres_  , i ) += band( repeat<R0,R1>( mat_    ), i );
               band( osres_ , i ) += band( repeat<R0,R1>( mat_    ), i );
               band( refres_, i ) += band( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( repeat<R0,R1>( omat_   ), i );
               band( odres_ , i ) += band( repeat<R0,R1>( omat_   ), i );
               band( sres_  , i ) += band( repeat<R0,R1>( omat_   ), i );
               band( osres_ , i ) += band( repeat<R0,R1>( omat_   ), i );
               band( refres_, i ) += band( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Band-wise repeat with addition assignment with evaluated matrix (runtime)
      {
         test_  = "Band-wise repeat with addition assignment with evaluated matrix (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( odres_ , i ) += band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( sres_  , i ) += band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( osres_ , i ) += band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( refres_, i ) += band( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( odres_ , i ) += band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( sres_  , i ) += band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( osres_ , i ) += band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( refres_, i ) += band( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Band-wise repeat with addition assignment with evaluated matrix (compile time)
      {
         test_  = "Band-wise repeat with addition assignment with evaluated matrix (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( odres_ , i ) += band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( sres_  , i ) += band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( osres_ , i ) += band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( refres_, i ) += band( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) += band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( odres_ , i ) += band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( sres_  , i ) += band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( osres_ , i ) += band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( refres_, i ) += band( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Band-wise repeat with subtraction assignment
      //=====================================================================================

      // Band-wise repeat with subtraction assignment with the given matrix (runtime)
      {
         test_  = "Band-wise repeat with subtraction assignment with the given matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( repeat( mat_   , R0, R1 ), i );
               band( odres_ , i ) -= band( repeat( mat_   , R0, R1 ), i );
               band( sres_  , i ) -= band( repeat( mat_   , R0, R1 ), i );
               band( osres_ , i ) -= band( repeat( mat_   , R0, R1 ), i );
               band( refres_, i ) -= band( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( repeat( omat_  , R0, R1 ), i );
               band( odres_ , i ) -= band( repeat( omat_  , R0, R1 ), i );
               band( sres_  , i ) -= band( repeat( omat_  , R0, R1 ), i );
               band( osres_ , i ) -= band( repeat( omat_  , R0, R1 ), i );
               band( refres_, i ) -= band( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Band-wise repeat with subtraction assignment with the given matrix (compile time)
      {
         test_  = "Band-wise repeat with subtraction assignment with the given matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( repeat<R0,R1>( mat_    ), i );
               band( odres_ , i ) -= band( repeat<R0,R1>( mat_    ), i );
               band( sres_  , i ) -= band( repeat<R0,R1>( mat_    ), i );
               band( osres_ , i ) -= band( repeat<R0,R1>( mat_    ), i );
               band( refres_, i ) -= band( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( repeat<R0,R1>( omat_   ), i );
               band( odres_ , i ) -= band( repeat<R0,R1>( omat_   ), i );
               band( sres_  , i ) -= band( repeat<R0,R1>( omat_   ), i );
               band( osres_ , i ) -= band( repeat<R0,R1>( omat_   ), i );
               band( refres_, i ) -= band( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Band-wise repeat with subtraction assignment with evaluated matrix (runtime)
      {
         test_  = "Band-wise repeat with subtraction assignment with evaluated matrix (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( odres_ , i ) -= band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( sres_  , i ) -= band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( osres_ , i ) -= band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( refres_, i ) -= band( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( odres_ , i ) -= band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( sres_  , i ) -= band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( osres_ , i ) -= band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( refres_, i ) -= band( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Band-wise repeat with subtraction assignment with evaluated matrix (compile time)
      {
         test_  = "Band-wise repeat with subtraction assignment with evaluated matrix (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( odres_ , i ) -= band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( sres_  , i ) -= band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( osres_ , i ) -= band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( refres_, i ) -= band( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) -= band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( odres_ , i ) -= band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( sres_  , i ) -= band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( osres_ , i ) -= band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( refres_, i ) -= band( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Band-wise repeat with multiplication assignment
      //=====================================================================================

      // Band-wise repeat with multiplication assignment with the given matrix (runtime)
      {
         test_  = "Band-wise repeat with multiplication assignment with the given matrix (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( repeat( mat_   , R0, R1 ), i );
               band( odres_ , i ) *= band( repeat( mat_   , R0, R1 ), i );
               band( sres_  , i ) *= band( repeat( mat_   , R0, R1 ), i );
               band( osres_ , i ) *= band( repeat( mat_   , R0, R1 ), i );
               band( refres_, i ) *= band( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( repeat( omat_  , R0, R1 ), i );
               band( odres_ , i ) *= band( repeat( omat_  , R0, R1 ), i );
               band( sres_  , i ) *= band( repeat( omat_  , R0, R1 ), i );
               band( osres_ , i ) *= band( repeat( omat_  , R0, R1 ), i );
               band( refres_, i ) *= band( repeat( refmat_, R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Band-wise repeat with multiplication assignment with the given matrix (compile time)
      {
         test_  = "Band-wise repeat with multiplication assignment with the given matrix (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( repeat<R0,R1>( mat_    ), i );
               band( odres_ , i ) *= band( repeat<R0,R1>( mat_    ), i );
               band( sres_  , i ) *= band( repeat<R0,R1>( mat_    ), i );
               band( osres_ , i ) *= band( repeat<R0,R1>( mat_    ), i );
               band( refres_, i ) *= band( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( repeat<R0,R1>( omat_   ), i );
               band( odres_ , i ) *= band( repeat<R0,R1>( omat_   ), i );
               band( sres_  , i ) *= band( repeat<R0,R1>( omat_   ), i );
               band( osres_ , i ) *= band( repeat<R0,R1>( omat_   ), i );
               band( refres_, i ) *= band( repeat<R0,R1>( refmat_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Band-wise repeat with multiplication assignment with evaluated matrix (runtime)
      {
         test_  = "Band-wise repeat with multiplication assignment with evaluated matrix (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( odres_ , i ) *= band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( sres_  , i ) *= band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( osres_ , i ) *= band( repeat( eval( mat_ )   , R0, R1 ), i );
               band( refres_, i ) *= band( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( odres_ , i ) *= band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( sres_  , i ) *= band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( osres_ , i ) *= band( repeat( eval( omat_ )  , R0, R1 ), i );
               band( refres_, i ) *= band( repeat( eval( refmat_ ), R0, R1 ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }

      // Band-wise repeat with multiplication assignment with evaluated matrix (compile time)
      {
         test_  = "Band-wise repeat with multiplication assignment with evaluated matrix (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( odres_ , i ) *= band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( sres_  , i ) *= band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( osres_ , i ) *= band( repeat<R0,R1>( eval( mat_ )    ), i );
               band( refres_, i ) *= band( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<MT>( ex );
         }

         checkResults<MT>();

         try {
            initResults();
            for( ptrdiff_t i=ibegin; i<iend; ++i ) {
               band( dres_  , i ) *= band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( odres_ , i ) *= band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( sres_  , i ) *= band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( osres_ , i ) *= band( repeat<R0,R1>( eval( omat_ )   ), i );
               band( refres_, i ) *= band( repeat<R0,R1>( eval( refmat_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<OMT>( ex );
         }

         checkResults<OMT>();
      }


      //=====================================================================================
      // Failure cases
      //=====================================================================================

      // Out-of-bounds access (invalid lower band index)
      {
         test_  = "Out-of-bounds band construction (invalid lower band index)";
         error_ = "Setup of out-of-bounds band succeeded";

         try {
            auto b = band( repeat( mat_, R0, R1 ), -mat_.rows()*R0 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( repeat<R0,R1>( mat_ ), -mat_.rows()*R0 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( repeat( omat_, R0, R1 ), -omat_.rows()*R0 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( repeat<R0,R1>( omat_ ), -omat_.rows()*R0 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
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
            auto b = band( repeat( mat_, R0, R1 ), mat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( repeat<R0,R1>( mat_ ), mat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( MT ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( repeat( omat_, R0, R1 ), omat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( repeat<R0,R1>( omat_ ), omat_.columns()*R1 );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Sparse matrix type:\n"
                << "     " << typeid( OMT ).name() << "\n"
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
/*!\brief Testing the customized sparse matrix repeat operation.
//
// \param op The custom operation to be tested.
// \param name The human-readable name of the operation.
// \return void
// \exception std::runtime_error Repeat error detected.
//
// This function tests the matrix repeat operation with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment in combination with
// a custom operation. In case any error resulting from the repeat operation or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT    // Type of the sparse matrix
        , size_t R0      // Compile time row-wise repetitions
        , size_t R1 >    // Compile time column-wise repetitions
template< typename OP >  // Type of the custom operation
void OperationTest<MT,R0,R1>::testCustomOperation( OP op, const std::string& name )
{
   using blaze::repeat;


   //=====================================================================================
   // Repeat operation
   //=====================================================================================

   // Customized repeat operation with the given matrix (runtime)
   {
      test_  = "Customized repeat operation with the given matrix (runtime)";
      error_ = "Failed repeat operation";

      try {
         initResults();
         dres_   = op( repeat( mat_, R0, R1 ) );
         odres_  = op( repeat( mat_, R0, R1 ) );
         sres_   = op( repeat( mat_, R0, R1 ) );
         osres_  = op( repeat( mat_, R0, R1 ) );
         refres_ = op( repeat( refmat_, R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   = op( repeat( omat_, R0, R1 ) );
         odres_  = op( repeat( omat_, R0, R1 ) );
         sres_   = op( repeat( omat_, R0, R1 ) );
         osres_  = op( repeat( omat_, R0, R1 ) );
         refres_ = op( repeat( refmat_, R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }

   // Customized repeat operation with the given matrix (compile time)
   {
      test_  = "Customized repeat operation with the given matrix (compile time)";
      error_ = "Failed repeat operation";

      try {
         initResults();
         dres_   = op( repeat<R0,R1>( mat_ ) );
         odres_  = op( repeat<R0,R1>( mat_ ) );
         sres_   = op( repeat<R0,R1>( mat_ ) );
         osres_  = op( repeat<R0,R1>( mat_ ) );
         refres_ = op( repeat<R0,R1>( refmat_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   = op( repeat<R0,R1>( omat_ ) );
         odres_  = op( repeat<R0,R1>( omat_ ) );
         sres_   = op( repeat<R0,R1>( omat_ ) );
         osres_  = op( repeat<R0,R1>( omat_ ) );
         refres_ = op( repeat<R0,R1>( refmat_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }

   // Customized repeat operation with evaluated matrix (runtime)
   {
      test_  = "Customized repeat operation with evaluated matrix (runtime)";
      error_ = "Failed repeat operation";

      try {
         initResults();
         dres_   = op( repeat( eval( mat_ ), R0, R1 ) );
         odres_  = op( repeat( eval( mat_ ), R0, R1 ) );
         sres_   = op( repeat( eval( mat_ ), R0, R1 ) );
         osres_  = op( repeat( eval( mat_ ), R0, R1 ) );
         refres_ = op( repeat( eval( refmat_ ), R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   = op( repeat( eval( omat_ ), R0, R1 ) );
         odres_  = op( repeat( eval( omat_ ), R0, R1 ) );
         sres_   = op( repeat( eval( omat_ ), R0, R1 ) );
         osres_  = op( repeat( eval( omat_ ), R0, R1 ) );
         refres_ = op( repeat( eval( refmat_ ), R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }

   // Customized repeat operation with evaluated matrix (compile time)
   {
      test_  = "Customized repeat operation with the given matrix (compile time)";
      error_ = "Failed repeat operation";

      try {
         initResults();
         dres_   = op( repeat<R0,R1>( eval( mat_ ) ) );
         odres_  = op( repeat<R0,R1>( eval( mat_ ) ) );
         sres_   = op( repeat<R0,R1>( eval( mat_ ) ) );
         osres_  = op( repeat<R0,R1>( eval( mat_ ) ) );
         refres_ = op( repeat<R0,R1>( eval( refmat_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   = op( repeat<R0,R1>( eval( omat_ ) ) );
         odres_  = op( repeat<R0,R1>( eval( omat_ ) ) );
         sres_   = op( repeat<R0,R1>( eval( omat_ ) ) );
         osres_  = op( repeat<R0,R1>( eval( omat_ ) ) );
         refres_ = op( repeat<R0,R1>( eval( refmat_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }


   //=====================================================================================
   // Repeat with addition assignment
   //=====================================================================================

   // Customized repeat with addition assignment with the given matrix (runtime)
   {
      test_  = "Customized repeat with addition assignment with the given matrix (runtime)";
      error_ = "Failed addition assignment";

      try {
         initResults();
         dres_   += op( repeat( mat_, R0, R1 ) );
         odres_  += op( repeat( mat_, R0, R1 ) );
         sres_   += op( repeat( mat_, R0, R1 ) );
         osres_  += op( repeat( mat_, R0, R1 ) );
         refres_ += op( repeat( refmat_, R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   += op( repeat( omat_, R0, R1 ) );
         odres_  += op( repeat( omat_, R0, R1 ) );
         sres_   += op( repeat( omat_, R0, R1 ) );
         osres_  += op( repeat( omat_, R0, R1 ) );
         refres_ += op( repeat( refmat_, R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }

   // Customized repeat with addition assignment with the given matrix (compile time)
   {
      test_  = "Customized repeat with addition assignment with the given matrix (compile time)";
      error_ = "Failed addition assignment";

      try {
         initResults();
         dres_   += op( repeat<R0,R1>( mat_ ) );
         odres_  += op( repeat<R0,R1>( mat_ ) );
         sres_   += op( repeat<R0,R1>( mat_ ) );
         osres_  += op( repeat<R0,R1>( mat_ ) );
         refres_ += op( repeat<R0,R1>( refmat_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   += op( repeat<R0,R1>( omat_ ) );
         odres_  += op( repeat<R0,R1>( omat_ ) );
         sres_   += op( repeat<R0,R1>( omat_ ) );
         osres_  += op( repeat<R0,R1>( omat_ ) );
         refres_ += op( repeat<R0,R1>( refmat_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }

   // Customized repeat with addition assignment with evaluated matrix (runtime)
   {
      test_  = "Customized repeat with addition assignment with evaluated matrix (runtime)";
      error_ = "Failed addition assignment";

      try {
         initResults();
         dres_   += op( repeat( eval( mat_ ), R0, R1 ) );
         odres_  += op( repeat( eval( mat_ ), R0, R1 ) );
         sres_   += op( repeat( eval( mat_ ), R0, R1 ) );
         osres_  += op( repeat( eval( mat_ ), R0, R1 ) );
         refres_ += op( repeat( eval( refmat_ ), R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   += op( repeat( eval( omat_ ), R0, R1 ) );
         odres_  += op( repeat( eval( omat_ ), R0, R1 ) );
         sres_   += op( repeat( eval( omat_ ), R0, R1 ) );
         osres_  += op( repeat( eval( omat_ ), R0, R1 ) );
         refres_ += op( repeat( eval( refmat_ ), R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }

   // Customized repeat with addition assignment with evaluated matrix (compile time)
   {
      test_  = "Customized repeat with addition assignment with the given matrix (compile time)";
      error_ = "Failed addition assignment";

      try {
         initResults();
         dres_   += op( repeat<R0,R1>( eval( mat_ ) ) );
         odres_  += op( repeat<R0,R1>( eval( mat_ ) ) );
         sres_   += op( repeat<R0,R1>( eval( mat_ ) ) );
         osres_  += op( repeat<R0,R1>( eval( mat_ ) ) );
         refres_ += op( repeat<R0,R1>( eval( refmat_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   += op( repeat<R0,R1>( eval( omat_ ) ) );
         odres_  += op( repeat<R0,R1>( eval( omat_ ) ) );
         sres_   += op( repeat<R0,R1>( eval( omat_ ) ) );
         osres_  += op( repeat<R0,R1>( eval( omat_ ) ) );
         refres_ += op( repeat<R0,R1>( eval( refmat_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }


   //=====================================================================================
   // Repeat with subtraction assignment
   //=====================================================================================

   // Customized repeat with subtraction assignment with the given matrix (runtime)
   {
      test_  = "Customized repeat with subtraction assignment with the given matrix (runtime)";
      error_ = "Failed subtraction assignment";

      try {
         initResults();
         dres_   -= op( repeat( mat_, R0, R1 ) );
         odres_  -= op( repeat( mat_, R0, R1 ) );
         sres_   -= op( repeat( mat_, R0, R1 ) );
         osres_  -= op( repeat( mat_, R0, R1 ) );
         refres_ -= op( repeat( refmat_, R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   -= op( repeat( omat_, R0, R1 ) );
         odres_  -= op( repeat( omat_, R0, R1 ) );
         sres_   -= op( repeat( omat_, R0, R1 ) );
         osres_  -= op( repeat( omat_, R0, R1 ) );
         refres_ -= op( repeat( refmat_, R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }

   // Customized repeat with subtraction assignment with the given matrix (compile time)
   {
      test_  = "Customized repeat with subtraction assignment with the given matrix (compile time)";
      error_ = "Failed subtraction assignment";

      try {
         initResults();
         dres_   -= op( repeat<R0,R1>( mat_ ) );
         odres_  -= op( repeat<R0,R1>( mat_ ) );
         sres_   -= op( repeat<R0,R1>( mat_ ) );
         osres_  -= op( repeat<R0,R1>( mat_ ) );
         refres_ -= op( repeat<R0,R1>( refmat_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   -= op( repeat<R0,R1>( omat_ ) );
         odres_  -= op( repeat<R0,R1>( omat_ ) );
         sres_   -= op( repeat<R0,R1>( omat_ ) );
         osres_  -= op( repeat<R0,R1>( omat_ ) );
         refres_ -= op( repeat<R0,R1>( refmat_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }

   // Customized repeat with subtraction assignment with evaluated matrix (runtime)
   {
      test_  = "Customized repeat with subtraction assignment with evaluated matrix (runtime)";
      error_ = "Failed subtraction assignment";

      try {
         initResults();
         dres_   -= op( repeat( eval( mat_ ), R0, R1 ) );
         odres_  -= op( repeat( eval( mat_ ), R0, R1 ) );
         sres_   -= op( repeat( eval( mat_ ), R0, R1 ) );
         osres_  -= op( repeat( eval( mat_ ), R0, R1 ) );
         refres_ -= op( repeat( eval( refmat_ ), R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   -= op( repeat( eval( omat_ ), R0, R1 ) );
         odres_  -= op( repeat( eval( omat_ ), R0, R1 ) );
         sres_   -= op( repeat( eval( omat_ ), R0, R1 ) );
         osres_  -= op( repeat( eval( omat_ ), R0, R1 ) );
         refres_ -= op( repeat( eval( refmat_ ), R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }

   // Customized repeat with subtraction assignment with evaluated matrix (compile time)
   {
      test_  = "Customized repeat with subtraction assignment with the given matrix (compile time)";
      error_ = "Failed subtraction assignment";

      try {
         initResults();
         dres_   -= op( repeat<R0,R1>( eval( mat_ ) ) );
         odres_  -= op( repeat<R0,R1>( eval( mat_ ) ) );
         sres_   -= op( repeat<R0,R1>( eval( mat_ ) ) );
         osres_  -= op( repeat<R0,R1>( eval( mat_ ) ) );
         refres_ -= op( repeat<R0,R1>( eval( refmat_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   -= op( repeat<R0,R1>( eval( omat_ ) ) );
         odres_  -= op( repeat<R0,R1>( eval( omat_ ) ) );
         sres_   -= op( repeat<R0,R1>( eval( omat_ ) ) );
         osres_  -= op( repeat<R0,R1>( eval( omat_ ) ) );
         refres_ -= op( repeat<R0,R1>( eval( refmat_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }


   //=====================================================================================
   // Repeat with Schur product assignment
   //=====================================================================================

   // Customized repeat with Schur product assignment with the given matrix (runtime)
   {
      test_  = "Customized repeat with Schur product assignment with the given matrix (runtime)";
      error_ = "Failed Schur product assignment";

      try {
         initResults();
         dres_   %= op( repeat( mat_, R0, R1 ) );
         odres_  %= op( repeat( mat_, R0, R1 ) );
         sres_   %= op( repeat( mat_, R0, R1 ) );
         osres_  %= op( repeat( mat_, R0, R1 ) );
         refres_ %= op( repeat( refmat_, R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   %= op( repeat( omat_, R0, R1 ) );
         odres_  %= op( repeat( omat_, R0, R1 ) );
         sres_   %= op( repeat( omat_, R0, R1 ) );
         osres_  %= op( repeat( omat_, R0, R1 ) );
         refres_ %= op( repeat( refmat_, R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }

   // Customized repeat with Schur product assignment with the given matrix (compile time)
   {
      test_  = "Customized repeat with Schur product assignment with the given matrix (compile time)";
      error_ = "Failed Schur product assignment";

      try {
         initResults();
         dres_   %= op( repeat<R0,R1>( mat_ ) );
         odres_  %= op( repeat<R0,R1>( mat_ ) );
         sres_   %= op( repeat<R0,R1>( mat_ ) );
         osres_  %= op( repeat<R0,R1>( mat_ ) );
         refres_ %= op( repeat<R0,R1>( refmat_ ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   %= op( repeat<R0,R1>( omat_ ) );
         odres_  %= op( repeat<R0,R1>( omat_ ) );
         sres_   %= op( repeat<R0,R1>( omat_ ) );
         osres_  %= op( repeat<R0,R1>( omat_ ) );
         refres_ %= op( repeat<R0,R1>( refmat_ ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }

   // Customized repeat with Schur product assignment with evaluated matrix (runtime)
   {
      test_  = "Customized repeat with Schur product assignment with evaluated matrix (runtime)";
      error_ = "Failed Schur product assignment";

      try {
         initResults();
         dres_   %= op( repeat( eval( mat_ ), R0, R1 ) );
         odres_  %= op( repeat( eval( mat_ ), R0, R1 ) );
         sres_   %= op( repeat( eval( mat_ ), R0, R1 ) );
         osres_  %= op( repeat( eval( mat_ ), R0, R1 ) );
         refres_ %= op( repeat( eval( refmat_ ), R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   %= op( repeat( eval( omat_ ), R0, R1 ) );
         odres_  %= op( repeat( eval( omat_ ), R0, R1 ) );
         sres_   %= op( repeat( eval( omat_ ), R0, R1 ) );
         osres_  %= op( repeat( eval( omat_ ), R0, R1 ) );
         refres_ %= op( repeat( eval( refmat_ ), R0, R1 ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
   }

   // Customized repeat with Schur product assignment with evaluated matrix (compile time)
   {
      test_  = "Customized repeat with Schur product assignment with the given matrix (compile time)";
      error_ = "Failed Schur product assignment";

      try {
         initResults();
         dres_   %= op( repeat<R0,R1>( eval( mat_ ) ) );
         odres_  %= op( repeat<R0,R1>( eval( mat_ ) ) );
         sres_   %= op( repeat<R0,R1>( eval( mat_ ) ) );
         osres_  %= op( repeat<R0,R1>( eval( mat_ ) ) );
         refres_ %= op( repeat<R0,R1>( eval( refmat_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<MT>( ex );
      }

      checkResults<MT>();

      try {
         initResults();
         dres_   %= op( repeat<R0,R1>( eval( omat_ ) ) );
         odres_  %= op( repeat<R0,R1>( eval( omat_ ) ) );
         sres_   %= op( repeat<R0,R1>( eval( omat_ ) ) );
         osres_  %= op( repeat<R0,R1>( eval( omat_ ) ) );
         refres_ %= op( repeat<R0,R1>( eval( refmat_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<OMT>( ex );
      }

      checkResults<OMT>();
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
template< typename MT   // Type of the sparse matrix
        , size_t R0     // Compile time row-wise repetitions
        , size_t R1 >   // Compile time column-wise repetitions
template< typename T >  // Type of the matrix operand
void OperationTest<MT,R0,R1>::checkResults()
{
   using blaze::IsRowMajorMatrix;

   if( !isEqual( dres_, refres_ ) || !isEqual( odres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result matrix detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse " << ( IsRowMajorMatrix<T>::value ? ( "row-major" ) : ( "column-major" ) ) << " matrix type:\n"
          << "     " << typeid( T ).name() << "\n"
          << "   Result:\n" << dres_ << "\n"
          << "   Result with opposite storage order:\n" << odres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( sres_, refres_ ) || !isEqual( osres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result matrix detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse " << ( IsRowMajorMatrix<T>::value ? ( "row-major" ) : ( "column-major" ) ) << " matrix type:\n"
          << "     " << typeid( T ).name() << "\n"
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
template< typename MT   // Type of the sparse matrix
        , size_t R0     // Compile time row-wise repetitions
        , size_t R1 >   // Compile time column-wise repetitions
template< typename T >  // Type of the matrix operand
void OperationTest<MT,R0,R1>::checkTransposeResults()
{
   using blaze::IsRowMajorMatrix;

   if( !isEqual( tdres_, refres_ ) || !isEqual( todres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result matrix detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse " << ( IsRowMajorMatrix<T>::value ? ( "row-major" ) : ( "column-major" ) ) << " matrix type:\n"
          << "     " << typeid( T ).name() << "\n"
          << "   Transpose result:\n" << tdres_ << "\n"
          << "   Transpose result with opposite storage order:\n" << todres_ << "\n"
          << "   Expected result:\n" << refres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( tsres_, refres_ ) || !isEqual( tosres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result matrix detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Sparse " << ( IsRowMajorMatrix<T>::value ? ( "row-major" ) : ( "column-major" ) ) << " matrix type:\n"
          << "     " << typeid( T ).name() << "\n"
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
template< typename MT   // Type of the sparse matrix
        , size_t R0     // Compile time row-wise repetitions
        , size_t R1 >   // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::checkExceptionMessage( const std::exception& ex, const std::string& message )
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
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::initResults()
{
   const blaze::UnderlyingBuiltin_t<DRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<DRE> max( randmax );

   resize( dres_, rows( mat_ )*R0, columns( mat_ )*R1 );
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
template< typename MT  // Type of the sparse matrix
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
void OperationTest<MT,R0,R1>::initTransposeResults()
{
   const blaze::UnderlyingBuiltin_t<TDRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<TDRE> max( randmax );

   resize( tdres_, columns( mat_ )*R1, rows( mat_ )*R0 );
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
template< typename MT   // Type of the sparse matrix
        , size_t R0     // Compile time row-wise repetitions
        , size_t R1 >   // Compile time column-wise repetitions
template< typename T >  // Type of the matrix operand
void OperationTest<MT,R0,R1>::convertException( const std::exception& ex )
{
   using blaze::IsRowMajorMatrix;

   std::ostringstream oss;
   oss << " Test : " << test_ << "\n"
       << " Error: " << error_ << "\n"
       << " Details:\n"
       << "   Random seed = " << blaze::getSeed() << "\n"
       << "   Sparse " << ( IsRowMajorMatrix<T>::value ? ( "row-major" ) : ( "column-major" ) ) << " matrix type:\n"
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
/*!\brief Testing the repeat operation for a specific matrix type.
//
// \param creator The creator for the sparse matrix.
// \return void
*/
template< typename MT >  // Type of the sparse matrix
void runTest( const Creator<MT>& creator )
{
   for( size_t rep=0UL; rep<BLAZETEST_REPETITIONS; ++rep ) {
      OperationTest<MT,3UL,9UL>{ creator };
      OperationTest<MT,6UL,6UL>{ creator };
      OperationTest<MT,9UL,3UL>{ creator };
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
/*!\brief Macro for the definition of a sparse matrix repeat operation test case.
*/
#define DEFINE_SMATREPEAT_OPERATION_TEST( MT ) \
   extern template class blazetest::mathtest::operations::smatrepeat::OperationTest<MT>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a sparse matrix repeat operation test case.
*/
#define RUN_SMATREPEAT_OPERATION_TEST( C ) \
   blazetest::mathtest::operations::smatrepeat::runTest( C )
/*! \endcond */
//*************************************************************************************************

} // namespace smatrepeat

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
