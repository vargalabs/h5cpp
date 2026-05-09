//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/dvecexpand/OperationTest.h
//  \brief Header file for the dense vector expansion operation test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_DVECEXPAND_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_DVECEXPAND_OPERATIONTEST_H_


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
#include <blaze/math/CompressedVector.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/traits/ExpandTrait.h>
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

namespace dvecexpand {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the dense vector expansion operation test.
//
// This class template represents one particular test of an expansion operation on a vector of
// a particular type. The template argument \a VT represents the type of the vector operand.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
class OperationTest
{
 private:
   //**Type definitions****************************************************************************
   using ET = blaze::ElementType_t<VT>;  //!< Element type.

   using TVT = blaze::TransposeType_t<VT>;  //!< Transpose vector type.

   using DRE   = blaze::ExpandTrait_t<VT,E>;    //!< Dense result type.
   using DET   = blaze::ElementType_t<DRE>;     //!< Element type of the dense result
   using ODRE  = blaze::OppositeType_t<DRE>;    //!< Dense result type with opposite storage order
   using TDRE  = blaze::TransposeType_t<DRE>;   //!< Transpose dense result type
   using TODRE = blaze::TransposeType_t<ODRE>;  //!< Transpose dense result type with opposite storage order

   using SRE   = blaze::CompressedMatrix<DET,true>;  //!< Sparse result type
   using SET   = blaze::ElementType_t<SRE>;          //!< Element type of the sparse result
   using OSRE  = blaze::OppositeType_t<SRE>;         //!< Sparse result type with opposite storage order
   using TSRE  = blaze::TransposeType_t<SRE>;        //!< Transpose sparse result type
   using TOSRE = blaze::TransposeType_t<OSRE>;       //!< Transpose sparse result type with opposite storage order

   using RT  = blaze::CompressedVector<ET,false>;  //!< Reference type.
   using RRE = blaze::ExpandTrait_t<RT,E>;         //!< Reference result type.

   using TRT  = blaze::TransposeType_t<RT>;   //!< Transpose reference type.
   using TRRE = blaze::ExpandTrait_t<TRT,E>;  //!< Transpose reference result type.

   //! Type of the vector expand expression (runtime argument)
   using VecExpandExprType1 =
      blaze::RemoveCVRef_t< decltype( blaze::expand( std::declval<VT>(), E ) ) >;

   //! Type of the vector expand expression (compile time argument)
   using VecExpandExprType2 =
      blaze::RemoveCVRef_t< decltype( blaze::expand<E>( std::declval<VT>() ) ) >;

   //! Type of the transpose vector expand expression (runtime argument)
   using TVecExpandExprType1 =
      blaze::RemoveCVRef_t< decltype( blaze::expand( std::declval<TVT>(), E ) ) >;

   //! Type of the transpose vector expand expression (compile time argument)
   using TVecExpandExprType2 =
      blaze::RemoveCVRef_t< decltype( blaze::expand<E>( std::declval<TVT>() ) ) >;
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
   VT    vec_;      //!< The dense vector operand.
   DRE   dres_;     //!< The dense result matrix.
   SRE   sres_;     //!< The sparse result matrix.
   ODRE  odres_;    //!< The dense result matrix with opposite storage order.
   OSRE  osres_;    //!< The sparse result matrix with opposite storage order.
   RT    refvec_;   //!< The reference vector.
   RRE   refres_;   //!< The reference result.
   TVT   tvec_;     //!< The transpose dense vector operand.
   TDRE  tdres_;    //!< The transpose dense result matrix.
   TSRE  tsres_;    //!< The transpose sparse result matrix.
   TODRE todres_;   //!< The transpose dense result matrix with opposite storage order.
   TOSRE tosres_;   //!< The transpose sparse result matrix with opposite storage order.
   TRT   trefvec_;  //!< The transpose reference vector.
   TRRE  trefres_;  //!< The transpose reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TVT   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( ODRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( TODRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OSRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TOSRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( RT    );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TRT   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RRE   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( TRRE  );

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE      ( VT    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE         ( TVT   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( DRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( ODRE  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TDRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TODRE );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( SRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( OSRE  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TSRE  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( TOSRE );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE      ( RT    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE         ( TRT   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( RRE   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( TRRE  );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DRE   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ODRE  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TDRE  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TODRE );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SRE   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OSRE  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TSRE  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TOSRE );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RRE   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TRRE  );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ET , blaze::ElementType_t<TVT>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<DRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<ODRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<TDRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<TODRE> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DET, blaze::ElementType_t<SRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<SRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<OSRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<TSRE>  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<TOSRE> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SET, blaze::ElementType_t<DRE>   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT , blaze::TransposeType_t<TVT> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RT , blaze::TransposeType_t<TRT> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( VecExpandExprType1, blaze::ResultType_t<VecExpandExprType1>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( VecExpandExprType1, blaze::TransposeType_t<VecExpandExprType1> );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( VecExpandExprType2, blaze::ResultType_t<VecExpandExprType2>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( VecExpandExprType2, blaze::TransposeType_t<VecExpandExprType2> );

   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( TVecExpandExprType1, blaze::ResultType_t<TVecExpandExprType1>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TVecExpandExprType1, blaze::TransposeType_t<TVecExpandExprType1> );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER     ( TVecExpandExprType2, blaze::ResultType_t<TVecExpandExprType2>    );
   BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_DIFFERENT_STORAGE_ORDER( TVecExpandExprType2, blaze::TransposeType_t<TVecExpandExprType2> );

   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( VecExpandExprType1 , blaze::BaseType_t<VecExpandExprType1>  );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( VecExpandExprType2 , blaze::BaseType_t<VecExpandExprType2>  );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TVecExpandExprType1, blaze::BaseType_t<TVecExpandExprType1> );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TVecExpandExprType2, blaze::BaseType_t<TVecExpandExprType2> );
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
/*!\brief Constructor for the dense vector expansion operation test.
//
// \param creator The creator for dense vector operand.
// \exception std::runtime_error Operation error detected.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
OperationTest<VT,E>::OperationTest( const Creator<VT>& creator )
   : vec_( creator() )     // The dense vector operand
   , dres_()               // The dense result matrix
   , sres_()               // The sparse result matrix
   , odres_()              // The dense result matrix with opposite storage order
   , osres_()              // The sparse result matrix with opposite storage order
   , refvec_( vec_ )       // The reference vector
   , refres_()             // The reference result
   , tvec_( trans(vec_) )  // The transpose dense vector operand
   , tdres_()              // The transpose dense result matrix
   , tsres_()              // The transpose sparse result matrix
   , todres_()             // The transpose dense result matrix with opposite storage order
   , tosres_()             // The transpose dense result matrix with opposite storage order
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
/*!\brief Tests on the initial status of the vector.
//
// \return void
// \exception std::runtime_error Initialization error detected.
//
// This function runs tests on the initial status of the vector. In case any initialization
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testInitialStatus()
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
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testAssignment()
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
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testEvaluation()
{
   using blaze::expand;


   //=====================================================================================
   // Testing the evaluation with a column vector
   //=====================================================================================

   {
      const auto res   ( evaluate( expand( vec_, E ) ) );
      const auto refres( evaluate( expand( refvec_, E ) ) );

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
      const auto res   ( evaluate( expand<E>( vec_ ) ) );
      const auto refres( evaluate( expand<E>( refvec_ ) ) );

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
      const auto res   ( evaluate( expand( eval( vec_ ), E ) ) );
      const auto refres( evaluate( expand( eval( refvec_ ), E ) ) );

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
      const auto res   ( evaluate( expand<E>( eval( vec_ ) ) ) );
      const auto refres( evaluate( expand<E>( eval( refvec_ ) ) ) );

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
      const auto res   ( evaluate( expand( tvec_, E ) ) );
      const auto refres( evaluate( expand( trefvec_, E ) ) );

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
      const auto res   ( evaluate( expand<E>( tvec_ ) ) );
      const auto refres( evaluate( expand<E>( trefvec_ ) ) );

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
      const auto res   ( evaluate( expand( eval( tvec_ ), E ) ) );
      const auto refres( evaluate( expand( eval( trefvec_ ), E ) ) );

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
      const auto res   ( evaluate( expand<E>( eval( tvec_ ) ) ) );
      const auto refres( evaluate( expand<E>( eval( trefvec_ ) ) ) );

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
/*!\brief Testing the matrix element access.
//
// \return void
// \exception std::runtime_error Element access error detected.
//
// This function tests the element access via the function call operator. In case any error
// is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testElementAccess()
{
   using blaze::equal;
   using blaze::expand;


   //=====================================================================================
   // Testing the element access with a column vector
   //=====================================================================================

   if( vec_.size() > 0UL && E > 0UL )
   {
      const size_t m( vec_.size() - 1UL );
      const size_t n( E-1UL );

      if( !equal( expand( vec_, E )(m,n), expand( refvec_, E )(m,n) ) ||
          !equal( expand( vec_, E ).at(m,n), expand( refvec_, E ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of expansion expression (runtime)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( expand<E>( vec_ )(m,n), expand<E>( refvec_ )(m,n) ) ||
          !equal( expand<E>( vec_ ).at(m,n), expand<E>( refvec_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of expansion expression (compile time)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( expand( eval( vec_ ), E )(m,n), expand( eval( refvec_ ), E )(m,n) ) ||
          !equal( expand( eval( vec_ ), E ).at(m,n), expand( eval( refvec_ ), E ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of evaluated expansion expression (runtime)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( expand<E>( eval( vec_ ) )(m,n), expand<E>( eval( refvec_ ) )(m,n) ) ||
          !equal( expand<E>( eval( vec_ ) ).at(m,n), expand<E>( eval( refvec_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of evaluated expansion expression (compile time)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense column vector type:\n"
             << "     " << typeid( VT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Testing the element access with a row vector
   //=====================================================================================

   if( tvec_.size() > 0UL && E > 0UL )
   {
      const size_t m( E-1UL );
      const size_t n( tvec_.size() - 1UL );

      if( !equal( expand( tvec_, E )(m,n), expand( trefvec_, E )(m,n) ) ||
          !equal( expand( tvec_, E ).at(m,n), expand( trefvec_, E ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of expansion expression (runtime)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense row vector type:\n"
             << "     " << typeid( TVT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( expand<E>( tvec_ )(m,n), expand<E>( trefvec_ )(m,n) ) ||
          !equal( expand<E>( tvec_ ).at(m,n), expand<E>( trefvec_ ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of expansion expression (compile time)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense row vector type:\n"
             << "     " << typeid( TVT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( expand( eval( tvec_ ), E )(m,n), expand( eval( trefvec_ ), E )(m,n) ) ||
          !equal( expand( eval( tvec_ ), E ).at(m,n), expand( eval( trefvec_ ), E ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of evaluated expansion expression (runtime)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense row vector type:\n"
             << "     " << typeid( TVT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( expand<E>( eval( tvec_ ) )(m,n), expand<E>( eval( trefvec_ ) )(m,n) ) ||
          !equal( expand<E>( eval( tvec_ ) ).at(m,n), expand<E>( eval( trefvec_ ) ).at(m,n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of evaluated expansion expression (compile time)\n"
             << " Error: Unequal resulting elements at element (" << m << "," << n << ") detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Dense row vector type:\n"
             << "     " << typeid( TVT ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the plain vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// expansion or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      using blaze::expand;


      //=====================================================================================
      // Expansion operation
      //=====================================================================================

      // Expansion operation with the given vector (runtime)
      {
         test_  = "Expansion operation with the given vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = expand( vec_, E );
            odres_  = expand( vec_, E );
            sres_   = expand( vec_, E );
            osres_  = expand( vec_, E );
            refres_ = expand( refvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = expand( tvec_, E );
            todres_  = expand( tvec_, E );
            tsres_   = expand( tvec_, E );
            tosres_  = expand( tvec_, E );
            trefres_ = expand( trefvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Expansion operation with the given vector (compile time)
      {
         test_  = "Expansion operation with the given vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = expand<E>( vec_ );
            odres_  = expand<E>( vec_ );
            sres_   = expand<E>( vec_ );
            osres_  = expand<E>( vec_ );
            refres_ = expand<E>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = expand<E>( tvec_ );
            todres_  = expand<E>( tvec_ );
            tsres_   = expand<E>( tvec_ );
            tosres_  = expand<E>( tvec_ );
            trefres_ = expand<E>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Expansion operation with evaluated vector (runtime)
      {
         test_  = "Expansion operation with evaluated vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = expand( eval( vec_ ), E );
            odres_  = expand( eval( vec_ ), E );
            sres_   = expand( eval( vec_ ), E );
            osres_  = expand( eval( vec_ ), E );
            refres_ = expand( eval( refvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = expand( eval( tvec_ ), E );
            todres_  = expand( eval( tvec_ ), E );
            tsres_   = expand( eval( tvec_ ), E );
            tosres_  = expand( eval( tvec_ ), E );
            trefres_ = expand( eval( trefvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Expansion operation with evaluated vector (compile time)
      {
         test_  = "Expansion operation with evaluated vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = expand<E>( eval( vec_ ) );
            odres_  = expand<E>( eval( vec_ ) );
            sres_   = expand<E>( eval( vec_ ) );
            osres_  = expand<E>( eval( vec_ ) );
            refres_ = expand<E>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = expand<E>( eval( tvec_ ) );
            todres_  = expand<E>( eval( tvec_ ) );
            tsres_   = expand<E>( eval( tvec_ ) );
            tosres_  = expand<E>( eval( tvec_ ) );
            trefres_ = expand<E>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Expansion with addition assignment
      //=====================================================================================

      // Expansion with addition assignment with the given vector (runtime)
      {
         test_  = "Expansion with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += expand( vec_, E );
            odres_  += expand( vec_, E );
            sres_   += expand( vec_, E );
            osres_  += expand( vec_, E );
            refres_ += expand( refvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += expand( tvec_, E );
            todres_  += expand( tvec_, E );
            tsres_   += expand( tvec_, E );
            tosres_  += expand( tvec_, E );
            trefres_ += expand( trefvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Expansion with addition assignment with the given vector (compile time)
      {
         test_  = "Expansion with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += expand<E>( vec_ );
            odres_  += expand<E>( vec_ );
            sres_   += expand<E>( vec_ );
            osres_  += expand<E>( vec_ );
            refres_ += expand<E>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += expand<E>( tvec_ );
            todres_  += expand<E>( tvec_ );
            tsres_   += expand<E>( tvec_ );
            tosres_  += expand<E>( tvec_ );
            trefres_ += expand<E>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Expansion with addition assignment with evaluated vector (runtime)
      {
         test_  = "Expansion with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += expand( eval( vec_ ), E );
            odres_  += expand( eval( vec_ ), E );
            sres_   += expand( eval( vec_ ), E );
            osres_  += expand( eval( vec_ ), E );
            refres_ += expand( eval( refvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += expand( eval( tvec_ ), E );
            todres_  += expand( eval( tvec_ ), E );
            tsres_   += expand( eval( tvec_ ), E );
            tosres_  += expand( eval( tvec_ ), E );
            trefres_ += expand( eval( trefvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Expansion with addition assignment with evaluated vector (compile time)
      {
         test_  = "Expansion with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += expand<E>( eval( vec_ ) );
            odres_  += expand<E>( eval( vec_ ) );
            sres_   += expand<E>( eval( vec_ ) );
            osres_  += expand<E>( eval( vec_ ) );
            refres_ += expand<E>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += expand<E>( eval( tvec_ ) );
            todres_  += expand<E>( eval( tvec_ ) );
            tsres_   += expand<E>( eval( tvec_ ) );
            tosres_  += expand<E>( eval( tvec_ ) );
            trefres_ += expand<E>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Expansion with subtraction assignment
      //=====================================================================================

      // Expansion with subtraction assignment with the given vector (runtime)
      {
         test_  = "Expansion with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= expand( vec_, E );
            odres_  -= expand( vec_, E );
            sres_   -= expand( vec_, E );
            osres_  -= expand( vec_, E );
            refres_ -= expand( refvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= expand( tvec_, E );
            todres_  -= expand( tvec_, E );
            tsres_   -= expand( tvec_, E );
            tosres_  -= expand( tvec_, E );
            trefres_ -= expand( trefvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Expansion with subtraction assignment with the given vector (compile time)
      {
         test_  = "Expansion with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= expand<E>( vec_ );
            odres_  -= expand<E>( vec_ );
            sres_   -= expand<E>( vec_ );
            osres_  -= expand<E>( vec_ );
            refres_ -= expand<E>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= expand<E>( tvec_ );
            todres_  -= expand<E>( tvec_ );
            tsres_   -= expand<E>( tvec_ );
            tosres_  -= expand<E>( tvec_ );
            trefres_ -= expand<E>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Expansion with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Expansion with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= expand( eval( vec_ ), E );
            odres_  -= expand( eval( vec_ ), E );
            sres_   -= expand( eval( vec_ ), E );
            osres_  -= expand( eval( vec_ ), E );
            refres_ -= expand( eval( refvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= expand( eval( tvec_ ), E );
            todres_  -= expand( eval( tvec_ ), E );
            tsres_   -= expand( eval( tvec_ ), E );
            tosres_  -= expand( eval( tvec_ ), E );
            trefres_ -= expand( eval( trefvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Expansion with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Expansion with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= expand<E>( eval( vec_ ) );
            odres_  -= expand<E>( eval( vec_ ) );
            sres_   -= expand<E>( eval( vec_ ) );
            osres_  -= expand<E>( eval( vec_ ) );
            refres_ -= expand<E>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= expand<E>( eval( tvec_ ) );
            todres_  -= expand<E>( eval( tvec_ ) );
            tsres_   -= expand<E>( eval( tvec_ ) );
            tosres_  -= expand<E>( eval( tvec_ ) );
            trefres_ -= expand<E>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Expansion with Schur product assignment
      //=====================================================================================

      // Expansion with Schur product assignment with the given vector (runtime)
      {
         test_  = "Expansion with Schur product assignment with the given vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= expand( vec_, E );
            odres_  %= expand( vec_, E );
            sres_   %= expand( vec_, E );
            osres_  %= expand( vec_, E );
            refres_ %= expand( refvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= expand( tvec_, E );
            todres_  %= expand( tvec_, E );
            tsres_   %= expand( tvec_, E );
            tosres_  %= expand( tvec_, E );
            trefres_ %= expand( trefvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Expansion with Schur product assignment with the given vector (compile time)
      {
         test_  = "Expansion with Schur product assignment with the given vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= expand<E>( vec_ );
            odres_  %= expand<E>( vec_ );
            sres_   %= expand<E>( vec_ );
            osres_  %= expand<E>( vec_ );
            refres_ %= expand<E>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= expand<E>( tvec_ );
            todres_  %= expand<E>( tvec_ );
            tsres_   %= expand<E>( tvec_ );
            tosres_  %= expand<E>( tvec_ );
            trefres_ %= expand<E>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Expansion with Schur product assignment with evaluated vector (runtime)
      {
         test_  = "Expansion with Schur product assignment with evaluated vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= expand( eval( vec_ ), E );
            odres_  %= expand( eval( vec_ ), E );
            sres_   %= expand( eval( vec_ ), E );
            osres_  %= expand( eval( vec_ ), E );
            refres_ %= expand( eval( refvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= expand( eval( tvec_ ), E );
            todres_  %= expand( eval( tvec_ ), E );
            tsres_   %= expand( eval( tvec_ ), E );
            tosres_  %= expand( eval( tvec_ ), E );
            trefres_ %= expand( eval( trefvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Expansion with Schur product assignment with evaluated vector (compile time)
      {
         test_  = "Expansion with Schur product assignment with evaluated vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= expand<E>( eval( vec_ ) );
            odres_  %= expand<E>( eval( vec_ ) );
            sres_   %= expand<E>( eval( vec_ ) );
            osres_  %= expand<E>( eval( vec_ ) );
            refres_ %= expand<E>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= expand<E>( eval( tvec_ ) );
            todres_  %= expand<E>( eval( tvec_ ) );
            tsres_   %= expand<E>( eval( tvec_ ) );
            tosres_  %= expand<E>( eval( tvec_ ) );
            trefres_ %= expand<E>( eval( trefvec_ ) );
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
/*!\brief Testing the negated dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the negated vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// expansion or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      using blaze::expand;


      //=====================================================================================
      // Negated expansion operation
      //=====================================================================================

      // Negated expansion operation with the given vector (runtime)
      {
         test_  = "Negated expansion operation with the given vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = -expand( vec_, E );
            odres_  = -expand( vec_, E );
            sres_   = -expand( vec_, E );
            osres_  = -expand( vec_, E );
            refres_ = -expand( refvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = -expand( tvec_, E );
            todres_  = -expand( tvec_, E );
            tsres_   = -expand( tvec_, E );
            tosres_  = -expand( tvec_, E );
            trefres_ = -expand( trefvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated expansion operation with the given vector (compile time)
      {
         test_  = "Negated expansion operation with the given vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = -expand<E>( vec_ );
            odres_  = -expand<E>( vec_ );
            sres_   = -expand<E>( vec_ );
            osres_  = -expand<E>( vec_ );
            refres_ = -expand<E>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = -expand<E>( tvec_ );
            todres_  = -expand<E>( tvec_ );
            tsres_   = -expand<E>( tvec_ );
            tosres_  = -expand<E>( tvec_ );
            trefres_ = -expand<E>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated expansion operation with evaluated vector (runtime)
      {
         test_  = "Negated expansion operation with evaluated vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = -expand( eval( vec_ ), E );
            odres_  = -expand( eval( vec_ ), E );
            sres_   = -expand( eval( vec_ ), E );
            osres_  = -expand( eval( vec_ ), E );
            refres_ = -expand( eval( refvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = -expand( eval( tvec_ ), E );
            todres_  = -expand( eval( tvec_ ), E );
            tsres_   = -expand( eval( tvec_ ), E );
            tosres_  = -expand( eval( tvec_ ), E );
            trefres_ = -expand( eval( trefvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated expansion operation with evaluated vector (compile time)
      {
         test_  = "Negated expansion operation with evaluated vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = -expand<E>( eval( vec_ ) );
            odres_  = -expand<E>( eval( vec_ ) );
            sres_   = -expand<E>( eval( vec_ ) );
            osres_  = -expand<E>( eval( vec_ ) );
            refres_ = -expand<E>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = -expand<E>( eval( tvec_ ) );
            todres_  = -expand<E>( eval( tvec_ ) );
            tsres_   = -expand<E>( eval( tvec_ ) );
            tosres_  = -expand<E>( eval( tvec_ ) );
            trefres_ = -expand<E>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Negated expansion with addition assignment
      //=====================================================================================

      // Negated expansion with addition assignment with the given vector (runtime)
      {
         test_  = "Negated expansion with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += -expand( vec_, E );
            odres_  += -expand( vec_, E );
            sres_   += -expand( vec_, E );
            osres_  += -expand( vec_, E );
            refres_ += -expand( refvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += -expand( tvec_, E );
            todres_  += -expand( tvec_, E );
            tsres_   += -expand( tvec_, E );
            tosres_  += -expand( tvec_, E );
            trefres_ += -expand( trefvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated expansion with addition assignment with the given vector (compile time)
      {
         test_  = "Negated expansion with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += -expand<E>( vec_ );
            odres_  += -expand<E>( vec_ );
            sres_   += -expand<E>( vec_ );
            osres_  += -expand<E>( vec_ );
            refres_ += -expand<E>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += -expand<E>( tvec_ );
            todres_  += -expand<E>( tvec_ );
            tsres_   += -expand<E>( tvec_ );
            tosres_  += -expand<E>( tvec_ );
            trefres_ += -expand<E>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated expansion with addition assignment with evaluated vector (runtime)
      {
         test_  = "Negated expansion with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += -expand( eval( vec_ ), E );
            odres_  += -expand( eval( vec_ ), E );
            sres_   += -expand( eval( vec_ ), E );
            osres_  += -expand( eval( vec_ ), E );
            refres_ += -expand( eval( refvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += -expand( eval( tvec_ ), E );
            todres_  += -expand( eval( tvec_ ), E );
            tsres_   += -expand( eval( tvec_ ), E );
            tosres_  += -expand( eval( tvec_ ), E );
            trefres_ += -expand( eval( trefvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated expansion with addition assignment with evaluated vector (compile time)
      {
         test_  = "Negated expansion with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += -expand<E>( eval( vec_ ) );
            odres_  += -expand<E>( eval( vec_ ) );
            sres_   += -expand<E>( eval( vec_ ) );
            osres_  += -expand<E>( eval( vec_ ) );
            refres_ += -expand<E>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += -expand<E>( eval( tvec_ ) );
            todres_  += -expand<E>( eval( tvec_ ) );
            tsres_   += -expand<E>( eval( tvec_ ) );
            tosres_  += -expand<E>( eval( tvec_ ) );
            trefres_ += -expand<E>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Negated expansion with subtraction assignment
      //=====================================================================================

      // Negated expansion with subtraction assignment with the given vector (runtime)
      {
         test_  = "Negated expansion with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= -expand( vec_, E );
            odres_  -= -expand( vec_, E );
            sres_   -= -expand( vec_, E );
            osres_  -= -expand( vec_, E );
            refres_ -= -expand( refvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= -expand( tvec_, E );
            todres_  -= -expand( tvec_, E );
            tsres_   -= -expand( tvec_, E );
            tosres_  -= -expand( tvec_, E );
            trefres_ -= -expand( trefvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated expansion with subtraction assignment with the given vector (compile time)
      {
         test_  = "Negated expansion with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= -expand<E>( vec_ );
            odres_  -= -expand<E>( vec_ );
            sres_   -= -expand<E>( vec_ );
            osres_  -= -expand<E>( vec_ );
            refres_ -= -expand<E>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= -expand<E>( tvec_ );
            todres_  -= -expand<E>( tvec_ );
            tsres_   -= -expand<E>( tvec_ );
            tosres_  -= -expand<E>( tvec_ );
            trefres_ -= -expand<E>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated expansion with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Negated expansion with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= -expand( eval( vec_ ), E );
            odres_  -= -expand( eval( vec_ ), E );
            sres_   -= -expand( eval( vec_ ), E );
            osres_  -= -expand( eval( vec_ ), E );
            refres_ -= -expand( eval( refvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= -expand( eval( tvec_ ), E );
            todres_  -= -expand( eval( tvec_ ), E );
            tsres_   -= -expand( eval( tvec_ ), E );
            tosres_  -= -expand( eval( tvec_ ), E );
            trefres_ -= -expand( eval( trefvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated expansion with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Negated expansion with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= -expand<E>( eval( vec_ ) );
            odres_  -= -expand<E>( eval( vec_ ) );
            sres_   -= -expand<E>( eval( vec_ ) );
            osres_  -= -expand<E>( eval( vec_ ) );
            refres_ -= -expand<E>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= -expand<E>( eval( tvec_ ) );
            todres_  -= -expand<E>( eval( tvec_ ) );
            tsres_   -= -expand<E>( eval( tvec_ ) );
            tosres_  -= -expand<E>( eval( tvec_ ) );
            trefres_ -= -expand<E>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Negated expansion with Schur product assignment
      //=====================================================================================

      // Negated expansion with Schur product assignment with the given vector (runtime)
      {
         test_  = "Negated expansion with Schur product assignment with the given vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= -expand( vec_, E );
            odres_  %= -expand( vec_, E );
            sres_   %= -expand( vec_, E );
            osres_  %= -expand( vec_, E );
            refres_ %= -expand( refvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= -expand( tvec_, E );
            todres_  %= -expand( tvec_, E );
            tsres_   %= -expand( tvec_, E );
            tosres_  %= -expand( tvec_, E );
            trefres_ %= -expand( trefvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated expansion with Schur product assignment with the given vector (compile time)
      {
         test_  = "Negated expansion with Schur product assignment with the given vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= -expand<E>( vec_ );
            odres_  %= -expand<E>( vec_ );
            sres_   %= -expand<E>( vec_ );
            osres_  %= -expand<E>( vec_ );
            refres_ %= -expand<E>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= -expand<E>( tvec_ );
            todres_  %= -expand<E>( tvec_ );
            tsres_   %= -expand<E>( tvec_ );
            tosres_  %= -expand<E>( tvec_ );
            trefres_ %= -expand<E>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated expansion with Schur product assignment with evaluated vector (runtime)
      {
         test_  = "Negated expansion with Schur product assignment with evaluated vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= -expand( eval( vec_ ), E );
            odres_  %= -expand( eval( vec_ ), E );
            sres_   %= -expand( eval( vec_ ), E );
            osres_  %= -expand( eval( vec_ ), E );
            refres_ %= -expand( eval( refvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= -expand( eval( tvec_ ), E );
            todres_  %= -expand( eval( tvec_ ), E );
            tsres_   %= -expand( eval( tvec_ ), E );
            tosres_  %= -expand( eval( tvec_ ), E );
            trefres_ %= -expand( eval( trefvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Negated expansion with Schur product assignment with evaluated vector (compile time)
      {
         test_  = "Negated expansion with Schur product assignment with evaluated vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= -expand<E>( eval( vec_ ) );
            odres_  %= -expand<E>( eval( vec_ ) );
            sres_   %= -expand<E>( eval( vec_ ) );
            osres_  %= -expand<E>( eval( vec_ ) );
            refres_ %= -expand<E>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= -expand<E>( eval( tvec_ ) );
            todres_  %= -expand<E>( eval( tvec_ ) );
            tsres_   %= -expand<E>( eval( tvec_ ) );
            tosres_  %= -expand<E>( eval( tvec_ ) );
            trefres_ %= -expand<E>( eval( trefvec_ ) );
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
/*!\brief Testing the scaled dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the scaled vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// expansion or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
template< typename T >    // Type of the scalar
void OperationTest<VT,E>::testScaledOperation( T scalar )
{
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE( T );

   if( scalar == T(0) )
      throw std::invalid_argument( "Invalid scalar parameter" );


#if BLAZETEST_MATHTEST_TEST_SCALED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SCALED_OPERATION > 1 )
   {
      using blaze::expand;


      //=====================================================================================
      // Scaled expansion (s*OP)
      //=====================================================================================

      // Scaled expansion operation with the given vector (s*OP, runtime)
      {
         test_  = "Scaled expansion operation with the given vector (s*OP, runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = scalar * expand( vec_, E );
            odres_  = scalar * expand( vec_, E );
            sres_   = scalar * expand( vec_, E );
            osres_  = scalar * expand( vec_, E );
            refres_ = scalar * expand( refvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = scalar * expand( tvec_, E );
            todres_  = scalar * expand( tvec_, E );
            tsres_   = scalar * expand( tvec_, E );
            tosres_  = scalar * expand( tvec_, E );
            trefres_ = scalar * expand( trefvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with the given vector (s*OP, compile time)
      {
         test_  = "Scaled expansion operation with the given vector (s*OP, compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = scalar * expand<E>( vec_ );
            odres_  = scalar * expand<E>( vec_ );
            sres_   = scalar * expand<E>( vec_ );
            osres_  = scalar * expand<E>( vec_ );
            refres_ = scalar * expand<E>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = scalar * expand<E>( tvec_ );
            todres_  = scalar * expand<E>( tvec_ );
            tsres_   = scalar * expand<E>( tvec_ );
            tosres_  = scalar * expand<E>( tvec_ );
            trefres_ = scalar * expand<E>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with evaluated vector (s*OP, runtime)
      {
         test_  = "Scaled expansion operation with evaluated vector (s*OP, runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = scalar * expand( eval( vec_ ), E );
            odres_  = scalar * expand( eval( vec_ ), E );
            sres_   = scalar * expand( eval( vec_ ), E );
            osres_  = scalar * expand( eval( vec_ ), E );
            refres_ = scalar * expand( eval( refvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = scalar * expand( eval( tvec_ ), E );
            todres_  = scalar * expand( eval( tvec_ ), E );
            tsres_   = scalar * expand( eval( tvec_ ), E );
            tosres_  = scalar * expand( eval( tvec_ ), E );
            trefres_ = scalar * expand( eval( trefvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with evaluated vector (s*OP, compile time)
      {
         test_  = "Scaled expansion operation with evaluated vector (s*OP, compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = scalar * expand<E>( eval( vec_ ) );
            odres_  = scalar * expand<E>( eval( vec_ ) );
            sres_   = scalar * expand<E>( eval( vec_ ) );
            osres_  = scalar * expand<E>( eval( vec_ ) );
            refres_ = scalar * expand<E>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = scalar * expand<E>( eval( tvec_ ) );
            todres_  = scalar * expand<E>( eval( tvec_ ) );
            tsres_   = scalar * expand<E>( eval( tvec_ ) );
            tosres_  = scalar * expand<E>( eval( tvec_ ) );
            trefres_ = scalar * expand<E>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled expansion (OP*s)
      //=====================================================================================

      // Scaled expansion operation with the given vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with the given vector (OP*s, runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = expand( vec_, E ) * scalar;
            odres_  = expand( vec_, E ) * scalar;
            sres_   = expand( vec_, E ) * scalar;
            osres_  = expand( vec_, E ) * scalar;
            refres_ = expand( refvec_, E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = expand( tvec_, E ) * scalar;
            todres_  = expand( tvec_, E ) * scalar;
            tsres_   = expand( tvec_, E ) * scalar;
            tosres_  = expand( tvec_, E ) * scalar;
            trefres_ = expand( trefvec_, E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with the given vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with the given vector (OP*s, compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = expand<E>( vec_ ) * scalar;
            odres_  = expand<E>( vec_ ) * scalar;
            sres_   = expand<E>( vec_ ) * scalar;
            osres_  = expand<E>( vec_ ) * scalar;
            refres_ = expand<E>( refvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = expand<E>( tvec_ ) * scalar;
            todres_  = expand<E>( tvec_ ) * scalar;
            tsres_   = expand<E>( tvec_ ) * scalar;
            tosres_  = expand<E>( tvec_ ) * scalar;
            trefres_ = expand<E>( trefvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with evaluated vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with evaluated vector (OP*s, runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = expand( eval( vec_ ), E ) * scalar;
            odres_  = expand( eval( vec_ ), E ) * scalar;
            sres_   = expand( eval( vec_ ), E ) * scalar;
            osres_  = expand( eval( vec_ ), E ) * scalar;
            refres_ = expand( eval( refvec_ ), E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = expand( eval( tvec_ ), E ) * scalar;
            todres_  = expand( eval( tvec_ ), E ) * scalar;
            tsres_   = expand( eval( tvec_ ), E ) * scalar;
            tosres_  = expand( eval( tvec_ ), E ) * scalar;
            trefres_ = expand( eval( trefvec_ ), E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with evaluated vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with evaluated vector (OP*s, compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = expand<E>( eval( vec_ ) ) * scalar;
            odres_  = expand<E>( eval( vec_ ) ) * scalar;
            sres_   = expand<E>( eval( vec_ ) ) * scalar;
            osres_  = expand<E>( eval( vec_ ) ) * scalar;
            refres_ = expand<E>( eval( refvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = expand<E>( eval( tvec_ ) ) * scalar;
            todres_  = expand<E>( eval( tvec_ ) ) * scalar;
            tsres_   = expand<E>( eval( tvec_ ) ) * scalar;
            tosres_  = expand<E>( eval( tvec_ ) ) * scalar;
            trefres_ = expand<E>( eval( trefvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled expansion (OP/s)
      //=====================================================================================

      // Scaled expansion operation with the given vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with the given vector (OP*s, runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = expand( vec_, E ) / scalar;
            odres_  = expand( vec_, E ) / scalar;
            sres_   = expand( vec_, E ) / scalar;
            osres_  = expand( vec_, E ) / scalar;
            refres_ = expand( refvec_, E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = expand( tvec_, E ) / scalar;
            todres_  = expand( tvec_, E ) / scalar;
            tsres_   = expand( tvec_, E ) / scalar;
            tosres_  = expand( tvec_, E ) / scalar;
            trefres_ = expand( trefvec_, E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with the given vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with the given vector (OP*s, compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = expand<E>( vec_ ) / scalar;
            odres_  = expand<E>( vec_ ) / scalar;
            sres_   = expand<E>( vec_ ) / scalar;
            osres_  = expand<E>( vec_ ) / scalar;
            refres_ = expand<E>( refvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = expand<E>( tvec_ ) / scalar;
            todres_  = expand<E>( tvec_ ) / scalar;
            tsres_   = expand<E>( tvec_ ) / scalar;
            tosres_  = expand<E>( tvec_ ) / scalar;
            trefres_ = expand<E>( trefvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with evaluated vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with evaluated vector (OP*s, runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = expand( eval( vec_ ), E ) / scalar;
            odres_  = expand( eval( vec_ ), E ) / scalar;
            sres_   = expand( eval( vec_ ), E ) / scalar;
            osres_  = expand( eval( vec_ ), E ) / scalar;
            refres_ = expand( eval( refvec_ ), E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = expand( eval( tvec_ ), E ) / scalar;
            todres_  = expand( eval( tvec_ ), E ) / scalar;
            tsres_   = expand( eval( tvec_ ), E ) / scalar;
            tosres_  = expand( eval( tvec_ ), E ) / scalar;
            trefres_ = expand( eval( trefvec_ ), E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with evaluated vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with evaluated vector (OP*s, compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            dres_   = expand<E>( eval( vec_ ) ) / scalar;
            odres_  = expand<E>( eval( vec_ ) ) / scalar;
            sres_   = expand<E>( eval( vec_ ) ) / scalar;
            osres_  = expand<E>( eval( vec_ ) ) / scalar;
            refres_ = expand<E>( eval( refvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   = expand<E>( eval( tvec_ ) ) / scalar;
            todres_  = expand<E>( eval( tvec_ ) ) / scalar;
            tsres_   = expand<E>( eval( tvec_ ) ) / scalar;
            tosres_  = expand<E>( eval( tvec_ ) ) / scalar;
            trefres_ = expand<E>( eval( trefvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled expansion with addition assignment (s*OP)
      //=====================================================================================

      // Scaled expansion operation with addition assignment with the given vector (s*OP, runtime)
      {
         test_  = "Scaled expansion operation with addition assignment with the given vector (s*OP, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += scalar * expand( vec_, E );
            odres_  += scalar * expand( vec_, E );
            sres_   += scalar * expand( vec_, E );
            osres_  += scalar * expand( vec_, E );
            refres_ += scalar * expand( refvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += scalar * expand( tvec_, E );
            todres_  += scalar * expand( tvec_, E );
            tsres_   += scalar * expand( tvec_, E );
            tosres_  += scalar * expand( tvec_, E );
            trefres_ += scalar * expand( trefvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with addition assignment with the given vector (s*OP, compile time)
      {
         test_  = "Scaled expansion operation with addition assignment with the given vector (s*OP, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += scalar * expand<E>( vec_ );
            odres_  += scalar * expand<E>( vec_ );
            sres_   += scalar * expand<E>( vec_ );
            osres_  += scalar * expand<E>( vec_ );
            refres_ += scalar * expand<E>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += scalar * expand<E>( tvec_ );
            todres_  += scalar * expand<E>( tvec_ );
            tsres_   += scalar * expand<E>( tvec_ );
            tosres_  += scalar * expand<E>( tvec_ );
            trefres_ += scalar * expand<E>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with addition assignment with evaluated vector (s*OP, runtime)
      {
         test_  = "Scaled expansion operation with addition assignment with evaluated vector (s*OP, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += scalar * expand( eval( vec_ ), E );
            odres_  += scalar * expand( eval( vec_ ), E );
            sres_   += scalar * expand( eval( vec_ ), E );
            osres_  += scalar * expand( eval( vec_ ), E );
            refres_ += scalar * expand( eval( refvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += scalar * expand( eval( tvec_ ), E );
            todres_  += scalar * expand( eval( tvec_ ), E );
            tsres_   += scalar * expand( eval( tvec_ ), E );
            tosres_  += scalar * expand( eval( tvec_ ), E );
            trefres_ += scalar * expand( eval( trefvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with addition assignment with evaluated vector (s*OP, compile time)
      {
         test_  = "Scaled expansion operation with addition assignment with evaluated vector (s*OP, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += scalar * expand<E>( eval( vec_ ) );
            odres_  += scalar * expand<E>( eval( vec_ ) );
            sres_   += scalar * expand<E>( eval( vec_ ) );
            osres_  += scalar * expand<E>( eval( vec_ ) );
            refres_ += scalar * expand<E>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += scalar * expand<E>( eval( tvec_ ) );
            todres_  += scalar * expand<E>( eval( tvec_ ) );
            tsres_   += scalar * expand<E>( eval( tvec_ ) );
            tosres_  += scalar * expand<E>( eval( tvec_ ) );
            trefres_ += scalar * expand<E>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled expansion with addition assignment (OP*s)
      //=====================================================================================

      // Scaled expansion operation with addition assignment with the given vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with addition assignment with the given vector (OP*s, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += expand( vec_, E ) * scalar;
            odres_  += expand( vec_, E ) * scalar;
            sres_   += expand( vec_, E ) * scalar;
            osres_  += expand( vec_, E ) * scalar;
            refres_ += expand( refvec_, E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += expand( tvec_, E ) * scalar;
            todres_  += expand( tvec_, E ) * scalar;
            tsres_   += expand( tvec_, E ) * scalar;
            tosres_  += expand( tvec_, E ) * scalar;
            trefres_ += expand( trefvec_, E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with addition assignment with the given vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with addition assignment with the given vector (OP*s, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += expand<E>( vec_ ) * scalar;
            odres_  += expand<E>( vec_ ) * scalar;
            sres_   += expand<E>( vec_ ) * scalar;
            osres_  += expand<E>( vec_ ) * scalar;
            refres_ += expand<E>( refvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += expand<E>( tvec_ ) * scalar;
            todres_  += expand<E>( tvec_ ) * scalar;
            tsres_   += expand<E>( tvec_ ) * scalar;
            tosres_  += expand<E>( tvec_ ) * scalar;
            trefres_ += expand<E>( trefvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with addition assignment with evaluated vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with addition assignment with evaluated vector (OP*s, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += expand( eval( vec_ ), E ) * scalar;
            odres_  += expand( eval( vec_ ), E ) * scalar;
            sres_   += expand( eval( vec_ ), E ) * scalar;
            osres_  += expand( eval( vec_ ), E ) * scalar;
            refres_ += expand( eval( refvec_ ), E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += expand( eval( tvec_ ), E ) * scalar;
            todres_  += expand( eval( tvec_ ), E ) * scalar;
            tsres_   += expand( eval( tvec_ ), E ) * scalar;
            tosres_  += expand( eval( tvec_ ), E ) * scalar;
            trefres_ += expand( eval( trefvec_ ), E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with addition assignment with evaluated vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with addition assignment with evaluated vector (OP*s, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += expand<E>( eval( vec_ ) ) * scalar;
            odres_  += expand<E>( eval( vec_ ) ) * scalar;
            sres_   += expand<E>( eval( vec_ ) ) * scalar;
            osres_  += expand<E>( eval( vec_ ) ) * scalar;
            refres_ += expand<E>( eval( refvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += expand<E>( eval( tvec_ ) ) * scalar;
            todres_  += expand<E>( eval( tvec_ ) ) * scalar;
            tsres_   += expand<E>( eval( tvec_ ) ) * scalar;
            tosres_  += expand<E>( eval( tvec_ ) ) * scalar;
            trefres_ += expand<E>( eval( trefvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled expansion (OP/s)
      //=====================================================================================

      // Scaled expansion operation with addition assignment with the given vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with addition assignment with the given vector (OP*s, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += expand( vec_, E ) / scalar;
            odres_  += expand( vec_, E ) / scalar;
            sres_   += expand( vec_, E ) / scalar;
            osres_  += expand( vec_, E ) / scalar;
            refres_ += expand( refvec_, E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += expand( tvec_, E ) / scalar;
            todres_  += expand( tvec_, E ) / scalar;
            tsres_   += expand( tvec_, E ) / scalar;
            tosres_  += expand( tvec_, E ) / scalar;
            trefres_ += expand( trefvec_, E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with addition assignment with the given vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with addition assignment with the given vector (OP*s, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += expand<E>( vec_ ) / scalar;
            odres_  += expand<E>( vec_ ) / scalar;
            sres_   += expand<E>( vec_ ) / scalar;
            osres_  += expand<E>( vec_ ) / scalar;
            refres_ += expand<E>( refvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += expand<E>( tvec_ ) / scalar;
            todres_  += expand<E>( tvec_ ) / scalar;
            tsres_   += expand<E>( tvec_ ) / scalar;
            tosres_  += expand<E>( tvec_ ) / scalar;
            trefres_ += expand<E>( trefvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with addition assignment with evaluated vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with addition assignment with evaluated vector (OP*s, runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += expand( eval( vec_ ), E ) / scalar;
            odres_  += expand( eval( vec_ ), E ) / scalar;
            sres_   += expand( eval( vec_ ), E ) / scalar;
            osres_  += expand( eval( vec_ ), E ) / scalar;
            refres_ += expand( eval( refvec_ ), E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += expand( eval( tvec_ ), E ) / scalar;
            todres_  += expand( eval( tvec_ ), E ) / scalar;
            tsres_   += expand( eval( tvec_ ), E ) / scalar;
            tosres_  += expand( eval( tvec_ ), E ) / scalar;
            trefres_ += expand( eval( trefvec_ ), E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with addition assignment with evaluated vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with addition assignment with evaluated vector (OP*s, compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            dres_   += expand<E>( eval( vec_ ) ) / scalar;
            odres_  += expand<E>( eval( vec_ ) ) / scalar;
            sres_   += expand<E>( eval( vec_ ) ) / scalar;
            osres_  += expand<E>( eval( vec_ ) ) / scalar;
            refres_ += expand<E>( eval( refvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   += expand<E>( eval( tvec_ ) ) / scalar;
            todres_  += expand<E>( eval( tvec_ ) ) / scalar;
            tsres_   += expand<E>( eval( tvec_ ) ) / scalar;
            tosres_  += expand<E>( eval( tvec_ ) ) / scalar;
            trefres_ += expand<E>( eval( trefvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled expansion with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled expansion operation with subtraction assignment with the given vector (s*OP, runtime)
      {
         test_  = "Scaled expansion operation with subtraction assignment with the given vector (s*OP, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= scalar * expand( vec_, E );
            odres_  -= scalar * expand( vec_, E );
            sres_   -= scalar * expand( vec_, E );
            osres_  -= scalar * expand( vec_, E );
            refres_ -= scalar * expand( refvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= scalar * expand( tvec_, E );
            todres_  -= scalar * expand( tvec_, E );
            tsres_   -= scalar * expand( tvec_, E );
            tosres_  -= scalar * expand( tvec_, E );
            trefres_ -= scalar * expand( trefvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with subtraction assignment with the given vector (s*OP, compile time)
      {
         test_  = "Scaled expansion operation with subtraction assignment with the given vector (s*OP, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= scalar * expand<E>( vec_ );
            odres_  -= scalar * expand<E>( vec_ );
            sres_   -= scalar * expand<E>( vec_ );
            osres_  -= scalar * expand<E>( vec_ );
            refres_ -= scalar * expand<E>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= scalar * expand<E>( tvec_ );
            todres_  -= scalar * expand<E>( tvec_ );
            tsres_   -= scalar * expand<E>( tvec_ );
            tosres_  -= scalar * expand<E>( tvec_ );
            trefres_ -= scalar * expand<E>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with subtraction assignment with evaluated vector (s*OP, runtime)
      {
         test_  = "Scaled expansion operation with subtraction assignment with evaluated vector (s*OP, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= scalar * expand( eval( vec_ ), E );
            odres_  -= scalar * expand( eval( vec_ ), E );
            sres_   -= scalar * expand( eval( vec_ ), E );
            osres_  -= scalar * expand( eval( vec_ ), E );
            refres_ -= scalar * expand( eval( refvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= scalar * expand( eval( tvec_ ), E );
            todres_  -= scalar * expand( eval( tvec_ ), E );
            tsres_   -= scalar * expand( eval( tvec_ ), E );
            tosres_  -= scalar * expand( eval( tvec_ ), E );
            trefres_ -= scalar * expand( eval( trefvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with subtraction assignment with evaluated vector (s*OP, compile time)
      {
         test_  = "Scaled expansion operation with subtraction assignment with evaluated vector (s*OP, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= scalar * expand<E>( eval( vec_ ) );
            odres_  -= scalar * expand<E>( eval( vec_ ) );
            sres_   -= scalar * expand<E>( eval( vec_ ) );
            osres_  -= scalar * expand<E>( eval( vec_ ) );
            refres_ -= scalar * expand<E>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= scalar * expand<E>( eval( tvec_ ) );
            todres_  -= scalar * expand<E>( eval( tvec_ ) );
            tsres_   -= scalar * expand<E>( eval( tvec_ ) );
            tosres_  -= scalar * expand<E>( eval( tvec_ ) );
            trefres_ -= scalar * expand<E>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled expansion with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled expansion operation with subtraction assignment with the given vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with subtraction assignment with the given vector (OP*s, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= expand( vec_, E ) * scalar;
            odres_  -= expand( vec_, E ) * scalar;
            sres_   -= expand( vec_, E ) * scalar;
            osres_  -= expand( vec_, E ) * scalar;
            refres_ -= expand( refvec_, E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= expand( tvec_, E ) * scalar;
            todres_  -= expand( tvec_, E ) * scalar;
            tsres_   -= expand( tvec_, E ) * scalar;
            tosres_  -= expand( tvec_, E ) * scalar;
            trefres_ -= expand( trefvec_, E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with subtraction assignment with the given vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with subtraction assignment with the given vector (OP*s, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= expand<E>( vec_ ) * scalar;
            odres_  -= expand<E>( vec_ ) * scalar;
            sres_   -= expand<E>( vec_ ) * scalar;
            osres_  -= expand<E>( vec_ ) * scalar;
            refres_ -= expand<E>( refvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= expand<E>( tvec_ ) * scalar;
            todres_  -= expand<E>( tvec_ ) * scalar;
            tsres_   -= expand<E>( tvec_ ) * scalar;
            tosres_  -= expand<E>( tvec_ ) * scalar;
            trefres_ -= expand<E>( trefvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with subtraction assignment with evaluated vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with subtraction assignment with evaluated vector (OP*s, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= expand( eval( vec_ ), E ) * scalar;
            odres_  -= expand( eval( vec_ ), E ) * scalar;
            sres_   -= expand( eval( vec_ ), E ) * scalar;
            osres_  -= expand( eval( vec_ ), E ) * scalar;
            refres_ -= expand( eval( refvec_ ), E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= expand( eval( tvec_ ), E ) * scalar;
            todres_  -= expand( eval( tvec_ ), E ) * scalar;
            tsres_   -= expand( eval( tvec_ ), E ) * scalar;
            tosres_  -= expand( eval( tvec_ ), E ) * scalar;
            trefres_ -= expand( eval( trefvec_ ), E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with subtraction assignment with evaluated vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with subtraction assignment with evaluated vector (OP*s, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= expand<E>( eval( vec_ ) ) * scalar;
            odres_  -= expand<E>( eval( vec_ ) ) * scalar;
            sres_   -= expand<E>( eval( vec_ ) ) * scalar;
            osres_  -= expand<E>( eval( vec_ ) ) * scalar;
            refres_ -= expand<E>( eval( refvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= expand<E>( eval( tvec_ ) ) * scalar;
            todres_  -= expand<E>( eval( tvec_ ) ) * scalar;
            tsres_   -= expand<E>( eval( tvec_ ) ) * scalar;
            tosres_  -= expand<E>( eval( tvec_ ) ) * scalar;
            trefres_ -= expand<E>( eval( trefvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled expansion (OP/s)
      //=====================================================================================

      // Scaled expansion operation with subtraction assignment with the given vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with subtraction assignment with the given vector (OP*s, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= expand( vec_, E ) / scalar;
            odres_  -= expand( vec_, E ) / scalar;
            sres_   -= expand( vec_, E ) / scalar;
            osres_  -= expand( vec_, E ) / scalar;
            refres_ -= expand( refvec_, E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= expand( tvec_, E ) / scalar;
            todres_  -= expand( tvec_, E ) / scalar;
            tsres_   -= expand( tvec_, E ) / scalar;
            tosres_  -= expand( tvec_, E ) / scalar;
            trefres_ -= expand( trefvec_, E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with subtraction assignment with the given vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with subtraction assignment with the given vector (OP*s, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= expand<E>( vec_ ) / scalar;
            odres_  -= expand<E>( vec_ ) / scalar;
            sres_   -= expand<E>( vec_ ) / scalar;
            osres_  -= expand<E>( vec_ ) / scalar;
            refres_ -= expand<E>( refvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= expand<E>( tvec_ ) / scalar;
            todres_  -= expand<E>( tvec_ ) / scalar;
            tsres_   -= expand<E>( tvec_ ) / scalar;
            tosres_  -= expand<E>( tvec_ ) / scalar;
            trefres_ -= expand<E>( trefvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with subtraction assignment with evaluated vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with subtraction assignment with evaluated vector (OP*s, runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= expand( eval( vec_ ), E ) / scalar;
            odres_  -= expand( eval( vec_ ), E ) / scalar;
            sres_   -= expand( eval( vec_ ), E ) / scalar;
            osres_  -= expand( eval( vec_ ), E ) / scalar;
            refres_ -= expand( eval( refvec_ ), E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= expand( eval( tvec_ ), E ) / scalar;
            todres_  -= expand( eval( tvec_ ), E ) / scalar;
            tsres_   -= expand( eval( tvec_ ), E ) / scalar;
            tosres_  -= expand( eval( tvec_ ), E ) / scalar;
            trefres_ -= expand( eval( trefvec_ ), E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with subtraction assignment with evaluated vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with subtraction assignment with evaluated vector (OP*s, compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            dres_   -= expand<E>( eval( vec_ ) ) / scalar;
            odres_  -= expand<E>( eval( vec_ ) ) / scalar;
            sres_   -= expand<E>( eval( vec_ ) ) / scalar;
            osres_  -= expand<E>( eval( vec_ ) ) / scalar;
            refres_ -= expand<E>( eval( refvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   -= expand<E>( eval( tvec_ ) ) / scalar;
            todres_  -= expand<E>( eval( tvec_ ) ) / scalar;
            tsres_   -= expand<E>( eval( tvec_ ) ) / scalar;
            tosres_  -= expand<E>( eval( tvec_ ) ) / scalar;
            trefres_ -= expand<E>( eval( trefvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled expansion with Schur product assignment (s*OP)
      //=====================================================================================

      // Scaled expansion operation with Schur product assignment with the given vector (s*OP, runtime)
      {
         test_  = "Scaled expansion operation with Schur product assignment with the given vector (s*OP, runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= scalar * expand( vec_, E );
            odres_  %= scalar * expand( vec_, E );
            sres_   %= scalar * expand( vec_, E );
            osres_  %= scalar * expand( vec_, E );
            refres_ %= scalar * expand( refvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= scalar * expand( tvec_, E );
            todres_  %= scalar * expand( tvec_, E );
            tsres_   %= scalar * expand( tvec_, E );
            tosres_  %= scalar * expand( tvec_, E );
            trefres_ %= scalar * expand( trefvec_, E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with Schur product assignment with the given vector (s*OP, compile time)
      {
         test_  = "Scaled expansion operation with Schur product assignment with the given vector (s*OP, compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= scalar * expand<E>( vec_ );
            odres_  %= scalar * expand<E>( vec_ );
            sres_   %= scalar * expand<E>( vec_ );
            osres_  %= scalar * expand<E>( vec_ );
            refres_ %= scalar * expand<E>( refvec_ );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= scalar * expand<E>( tvec_ );
            todres_  %= scalar * expand<E>( tvec_ );
            tsres_   %= scalar * expand<E>( tvec_ );
            tosres_  %= scalar * expand<E>( tvec_ );
            trefres_ %= scalar * expand<E>( trefvec_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with Schur product assignment with evaluated vector (s*OP, runtime)
      {
         test_  = "Scaled expansion operation with Schur product assignment with evaluated vector (s*OP, runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= scalar * expand( eval( vec_ ), E );
            odres_  %= scalar * expand( eval( vec_ ), E );
            sres_   %= scalar * expand( eval( vec_ ), E );
            osres_  %= scalar * expand( eval( vec_ ), E );
            refres_ %= scalar * expand( eval( refvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= scalar * expand( eval( tvec_ ), E );
            todres_  %= scalar * expand( eval( tvec_ ), E );
            tsres_   %= scalar * expand( eval( tvec_ ), E );
            tosres_  %= scalar * expand( eval( tvec_ ), E );
            trefres_ %= scalar * expand( eval( trefvec_ ), E );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with Schur product assignment with evaluated vector (s*OP, compile time)
      {
         test_  = "Scaled expansion operation with Schur product assignment with evaluated vector (s*OP, compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= scalar * expand<E>( eval( vec_ ) );
            odres_  %= scalar * expand<E>( eval( vec_ ) );
            sres_   %= scalar * expand<E>( eval( vec_ ) );
            osres_  %= scalar * expand<E>( eval( vec_ ) );
            refres_ %= scalar * expand<E>( eval( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= scalar * expand<E>( eval( tvec_ ) );
            todres_  %= scalar * expand<E>( eval( tvec_ ) );
            tsres_   %= scalar * expand<E>( eval( tvec_ ) );
            tosres_  %= scalar * expand<E>( eval( tvec_ ) );
            trefres_ %= scalar * expand<E>( eval( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled expansion with Schur product assignment (OP*s)
      //=====================================================================================

      // Scaled expansion operation with Schur product assignment with the given vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with Schur product assignment with the given vector (OP*s, runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= expand( vec_, E ) * scalar;
            odres_  %= expand( vec_, E ) * scalar;
            sres_   %= expand( vec_, E ) * scalar;
            osres_  %= expand( vec_, E ) * scalar;
            refres_ %= expand( refvec_, E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= expand( tvec_, E ) * scalar;
            todres_  %= expand( tvec_, E ) * scalar;
            tsres_   %= expand( tvec_, E ) * scalar;
            tosres_  %= expand( tvec_, E ) * scalar;
            trefres_ %= expand( trefvec_, E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with Schur product assignment with the given vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with Schur product assignment with the given vector (OP*s, compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= expand<E>( vec_ ) * scalar;
            odres_  %= expand<E>( vec_ ) * scalar;
            sres_   %= expand<E>( vec_ ) * scalar;
            osres_  %= expand<E>( vec_ ) * scalar;
            refres_ %= expand<E>( refvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= expand<E>( tvec_ ) * scalar;
            todres_  %= expand<E>( tvec_ ) * scalar;
            tsres_   %= expand<E>( tvec_ ) * scalar;
            tosres_  %= expand<E>( tvec_ ) * scalar;
            trefres_ %= expand<E>( trefvec_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with Schur product assignment with evaluated vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with Schur product assignment with evaluated vector (OP*s, runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= expand( eval( vec_ ), E ) * scalar;
            odres_  %= expand( eval( vec_ ), E ) * scalar;
            sres_   %= expand( eval( vec_ ), E ) * scalar;
            osres_  %= expand( eval( vec_ ), E ) * scalar;
            refres_ %= expand( eval( refvec_ ), E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= expand( eval( tvec_ ), E ) * scalar;
            todres_  %= expand( eval( tvec_ ), E ) * scalar;
            tsres_   %= expand( eval( tvec_ ), E ) * scalar;
            tosres_  %= expand( eval( tvec_ ), E ) * scalar;
            trefres_ %= expand( eval( trefvec_ ), E ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with Schur product assignment with evaluated vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with Schur product assignment with evaluated vector (OP*s, compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= expand<E>( eval( vec_ ) ) * scalar;
            odres_  %= expand<E>( eval( vec_ ) ) * scalar;
            sres_   %= expand<E>( eval( vec_ ) ) * scalar;
            osres_  %= expand<E>( eval( vec_ ) ) * scalar;
            refres_ %= expand<E>( eval( refvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= expand<E>( eval( tvec_ ) ) * scalar;
            todres_  %= expand<E>( eval( tvec_ ) ) * scalar;
            tsres_   %= expand<E>( eval( tvec_ ) ) * scalar;
            tosres_  %= expand<E>( eval( tvec_ ) ) * scalar;
            trefres_ %= expand<E>( eval( trefvec_ ) ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Scaled expansion (OP/s)
      //=====================================================================================

      // Scaled expansion operation with Schur product assignment with the given vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with Schur product assignment with the given vector (OP*s, runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= expand( vec_, E ) / scalar;
            odres_  %= expand( vec_, E ) / scalar;
            sres_   %= expand( vec_, E ) / scalar;
            osres_  %= expand( vec_, E ) / scalar;
            refres_ %= expand( refvec_, E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= expand( tvec_, E ) / scalar;
            todres_  %= expand( tvec_, E ) / scalar;
            tsres_   %= expand( tvec_, E ) / scalar;
            tosres_  %= expand( tvec_, E ) / scalar;
            trefres_ %= expand( trefvec_, E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with Schur product assignment with the given vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with Schur product assignment with the given vector (OP*s, compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= expand<E>( vec_ ) / scalar;
            odres_  %= expand<E>( vec_ ) / scalar;
            sres_   %= expand<E>( vec_ ) / scalar;
            osres_  %= expand<E>( vec_ ) / scalar;
            refres_ %= expand<E>( refvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= expand<E>( tvec_ ) / scalar;
            todres_  %= expand<E>( tvec_ ) / scalar;
            tsres_   %= expand<E>( tvec_ ) / scalar;
            tosres_  %= expand<E>( tvec_ ) / scalar;
            trefres_ %= expand<E>( trefvec_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with Schur product assignment with evaluated vector (OP*s, runtime)
      {
         test_  = "Scaled expansion operation with Schur product assignment with evaluated vector (OP*s, runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= expand( eval( vec_ ), E ) / scalar;
            odres_  %= expand( eval( vec_ ), E ) / scalar;
            sres_   %= expand( eval( vec_ ), E ) / scalar;
            osres_  %= expand( eval( vec_ ), E ) / scalar;
            refres_ %= expand( eval( refvec_ ), E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= expand( eval( tvec_ ), E ) / scalar;
            todres_  %= expand( eval( tvec_ ), E ) / scalar;
            tsres_   %= expand( eval( tvec_ ), E ) / scalar;
            tosres_  %= expand( eval( tvec_ ), E ) / scalar;
            trefres_ %= expand( eval( trefvec_ ), E ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Scaled expansion operation with Schur product assignment with evaluated vector (OP*s, compile time)
      {
         test_  = "Scaled expansion operation with Schur product assignment with evaluated vector (OP*s, compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            dres_   %= expand<E>( eval( vec_ ) ) / scalar;
            odres_  %= expand<E>( eval( vec_ ) ) / scalar;
            sres_   %= expand<E>( eval( vec_ ) ) / scalar;
            osres_  %= expand<E>( eval( vec_ ) ) / scalar;
            refres_ %= expand<E>( eval( refvec_ ) ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            tdres_   %= expand<E>( eval( tvec_ ) ) / scalar;
            todres_  %= expand<E>( eval( tvec_ ) ) / scalar;
            tsres_   %= expand<E>( eval( tvec_ ) ) / scalar;
            tosres_  %= expand<E>( eval( tvec_ ) ) / scalar;
            trefres_ %= expand<E>( eval( trefvec_ ) ) / scalar;
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
/*!\brief Testing the transpose dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the transpose vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// expansion or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_TRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_TRANS_OPERATION > 1 )
   {
      using blaze::expand;


      //=====================================================================================
      // Transpose expansion operation
      //=====================================================================================

      // Transpose expansion operation with the given vector (runtime)
      {
         test_  = "Transpose expansion operation with the given vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initTransposeResults();
            tdres_   = trans( expand( vec_, E ) );
            todres_  = trans( expand( vec_, E ) );
            tsres_   = trans( expand( vec_, E ) );
            tosres_  = trans( expand( vec_, E ) );
            trefres_ = trans( expand( refvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = trans( expand( tvec_, E ) );
            odres_  = trans( expand( tvec_, E ) );
            sres_   = trans( expand( tvec_, E ) );
            osres_  = trans( expand( tvec_, E ) );
            refres_ = trans( expand( trefvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose expansion operation with the given vector (compile time)
      {
         test_  = "Transpose expansion operation with the given vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initTransposeResults();
            tdres_   = trans( expand<E>( vec_ ) );
            todres_  = trans( expand<E>( vec_ ) );
            tsres_   = trans( expand<E>( vec_ ) );
            tosres_  = trans( expand<E>( vec_ ) );
            trefres_ = trans( expand<E>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = trans( expand<E>( tvec_ ) );
            odres_  = trans( expand<E>( tvec_ ) );
            sres_   = trans( expand<E>( tvec_ ) );
            osres_  = trans( expand<E>( tvec_ ) );
            refres_ = trans( expand<E>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose expansion operation with evaluated vector (runtime)
      {
         test_  = "Transpose expansion operation with evaluated vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initTransposeResults();
            tdres_   = trans( expand( eval( vec_ ), E ) );
            todres_  = trans( expand( eval( vec_ ), E ) );
            tsres_   = trans( expand( eval( vec_ ), E ) );
            tosres_  = trans( expand( eval( vec_ ), E ) );
            trefres_ = trans( expand( eval( refvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = trans( expand( eval( tvec_ ), E ) );
            odres_  = trans( expand( eval( tvec_ ), E ) );
            sres_   = trans( expand( eval( tvec_ ), E ) );
            osres_  = trans( expand( eval( tvec_ ), E ) );
            refres_ = trans( expand( eval( trefvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose expansion operation with evaluated vector (compile time)
      {
         test_  = "Transpose expansion operation with evaluated vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initTransposeResults();
            tdres_   = trans( expand<E>( eval( vec_ ) ) );
            todres_  = trans( expand<E>( eval( vec_ ) ) );
            tsres_   = trans( expand<E>( eval( vec_ ) ) );
            tosres_  = trans( expand<E>( eval( vec_ ) ) );
            trefres_ = trans( expand<E>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = trans( expand<E>( eval( tvec_ ) ) );
            odres_  = trans( expand<E>( eval( tvec_ ) ) );
            sres_   = trans( expand<E>( eval( tvec_ ) ) );
            osres_  = trans( expand<E>( eval( tvec_ ) ) );
            refres_ = trans( expand<E>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }


      //=====================================================================================
      // Transpose expansion with addition assignment
      //=====================================================================================

      // Transpose expansion with addition assignment with the given vector (runtime)
      {
         test_  = "Transpose expansion with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += trans( expand( vec_, E ) );
            todres_  += trans( expand( vec_, E ) );
            tsres_   += trans( expand( vec_, E ) );
            tosres_  += trans( expand( vec_, E ) );
            trefres_ += trans( expand( refvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += trans( expand( tvec_, E ) );
            odres_  += trans( expand( tvec_, E ) );
            sres_   += trans( expand( tvec_, E ) );
            osres_  += trans( expand( tvec_, E ) );
            refres_ += trans( expand( trefvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose expansion with addition assignment with the given vector (compile time)
      {
         test_  = "Transpose expansion with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += trans( expand<E>( vec_ ) );
            todres_  += trans( expand<E>( vec_ ) );
            tsres_   += trans( expand<E>( vec_ ) );
            tosres_  += trans( expand<E>( vec_ ) );
            trefres_ += trans( expand<E>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += trans( expand<E>( tvec_ ) );
            odres_  += trans( expand<E>( tvec_ ) );
            sres_   += trans( expand<E>( tvec_ ) );
            osres_  += trans( expand<E>( tvec_ ) );
            refres_ += trans( expand<E>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose expansion with addition assignment with evaluated vector (runtime)
      {
         test_  = "Transpose expansion with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += trans( expand( eval( vec_ ), E ) );
            todres_  += trans( expand( eval( vec_ ), E ) );
            tsres_   += trans( expand( eval( vec_ ), E ) );
            tosres_  += trans( expand( eval( vec_ ), E ) );
            trefres_ += trans( expand( eval( refvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += trans( expand( eval( tvec_ ), E ) );
            odres_  += trans( expand( eval( tvec_ ), E ) );
            sres_   += trans( expand( eval( tvec_ ), E ) );
            osres_  += trans( expand( eval( tvec_ ), E ) );
            refres_ += trans( expand( eval( trefvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose expansion with addition assignment with evaluated vector (compile time)
      {
         test_  = "Transpose expansion with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += trans( expand<E>( eval( vec_ ) ) );
            todres_  += trans( expand<E>( eval( vec_ ) ) );
            tsres_   += trans( expand<E>( eval( vec_ ) ) );
            tosres_  += trans( expand<E>( eval( vec_ ) ) );
            trefres_ += trans( expand<E>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += trans( expand<E>( eval( tvec_ ) ) );
            odres_  += trans( expand<E>( eval( tvec_ ) ) );
            sres_   += trans( expand<E>( eval( tvec_ ) ) );
            osres_  += trans( expand<E>( eval( tvec_ ) ) );
            refres_ += trans( expand<E>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }


      //=====================================================================================
      // Transpose expansion with subtraction assignment
      //=====================================================================================

      // Transpose expansion with subtraction assignment with the given vector (runtime)
      {
         test_  = "Transpose expansion with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= trans( expand( vec_, E ) );
            todres_  -= trans( expand( vec_, E ) );
            tsres_   -= trans( expand( vec_, E ) );
            tosres_  -= trans( expand( vec_, E ) );
            trefres_ -= trans( expand( refvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= trans( expand( tvec_, E ) );
            odres_  -= trans( expand( tvec_, E ) );
            sres_   -= trans( expand( tvec_, E ) );
            osres_  -= trans( expand( tvec_, E ) );
            refres_ -= trans( expand( trefvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose expansion with subtraction assignment with the given vector (compile time)
      {
         test_  = "Transpose expansion with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= trans( expand<E>( vec_ ) );
            todres_  -= trans( expand<E>( vec_ ) );
            tsres_   -= trans( expand<E>( vec_ ) );
            tosres_  -= trans( expand<E>( vec_ ) );
            trefres_ -= trans( expand<E>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= trans( expand<E>( tvec_ ) );
            odres_  -= trans( expand<E>( tvec_ ) );
            sres_   -= trans( expand<E>( tvec_ ) );
            osres_  -= trans( expand<E>( tvec_ ) );
            refres_ -= trans( expand<E>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose expansion with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Transpose expansion with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= trans( expand( eval( vec_ ), E ) );
            todres_  -= trans( expand( eval( vec_ ), E ) );
            tsres_   -= trans( expand( eval( vec_ ), E ) );
            tosres_  -= trans( expand( eval( vec_ ), E ) );
            trefres_ -= trans( expand( eval( refvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= trans( expand( eval( tvec_ ), E ) );
            odres_  -= trans( expand( eval( tvec_ ), E ) );
            sres_   -= trans( expand( eval( tvec_ ), E ) );
            osres_  -= trans( expand( eval( tvec_ ), E ) );
            refres_ -= trans( expand( eval( trefvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose expansion with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Transpose expansion with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= trans( expand<E>( eval( vec_ ) ) );
            todres_  -= trans( expand<E>( eval( vec_ ) ) );
            tsres_   -= trans( expand<E>( eval( vec_ ) ) );
            tosres_  -= trans( expand<E>( eval( vec_ ) ) );
            trefres_ -= trans( expand<E>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= trans( expand<E>( eval( tvec_ ) ) );
            odres_  -= trans( expand<E>( eval( tvec_ ) ) );
            sres_   -= trans( expand<E>( eval( tvec_ ) ) );
            osres_  -= trans( expand<E>( eval( tvec_ ) ) );
            refres_ -= trans( expand<E>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }


      //=====================================================================================
      // Transpose expansion with Schur product assignment
      //=====================================================================================

      // Transpose expansion with Schur product assignment with the given vector (runtime)
      {
         test_  = "Transpose expansion with Schur product assignment with the given vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initTransposeResults();
            tdres_   %= trans( expand( vec_, E ) );
            todres_  %= trans( expand( vec_, E ) );
            tsres_   %= trans( expand( vec_, E ) );
            tosres_  %= trans( expand( vec_, E ) );
            trefres_ %= trans( expand( refvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   %= trans( expand( tvec_, E ) );
            odres_  %= trans( expand( tvec_, E ) );
            sres_   %= trans( expand( tvec_, E ) );
            osres_  %= trans( expand( tvec_, E ) );
            refres_ %= trans( expand( trefvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose expansion with Schur product assignment with the given vector (compile time)
      {
         test_  = "Transpose expansion with Schur product assignment with the given vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initTransposeResults();
            tdres_   %= trans( expand<E>( vec_ ) );
            todres_  %= trans( expand<E>( vec_ ) );
            tsres_   %= trans( expand<E>( vec_ ) );
            tosres_  %= trans( expand<E>( vec_ ) );
            trefres_ %= trans( expand<E>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   %= trans( expand<E>( tvec_ ) );
            odres_  %= trans( expand<E>( tvec_ ) );
            sres_   %= trans( expand<E>( tvec_ ) );
            osres_  %= trans( expand<E>( tvec_ ) );
            refres_ %= trans( expand<E>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose expansion with Schur product assignment with evaluated vector (runtime)
      {
         test_  = "Transpose expansion with Schur product assignment with evaluated vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initTransposeResults();
            tdres_   %= trans( expand( eval( vec_ ), E ) );
            todres_  %= trans( expand( eval( vec_ ), E ) );
            tsres_   %= trans( expand( eval( vec_ ), E ) );
            tosres_  %= trans( expand( eval( vec_ ), E ) );
            trefres_ %= trans( expand( eval( refvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   %= trans( expand( eval( tvec_ ), E ) );
            odres_  %= trans( expand( eval( tvec_ ), E ) );
            sres_   %= trans( expand( eval( tvec_ ), E ) );
            osres_  %= trans( expand( eval( tvec_ ), E ) );
            refres_ %= trans( expand( eval( trefvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Transpose expansion with Schur product assignment with evaluated vector (compile time)
      {
         test_  = "Transpose expansion with Schur product assignment with evaluated vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initTransposeResults();
            tdres_   %= trans( expand<E>( eval( vec_ ) ) );
            todres_  %= trans( expand<E>( eval( vec_ ) ) );
            tsres_   %= trans( expand<E>( eval( vec_ ) ) );
            tosres_  %= trans( expand<E>( eval( vec_ ) ) );
            trefres_ %= trans( expand<E>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   %= trans( expand<E>( eval( tvec_ ) ) );
            odres_  %= trans( expand<E>( eval( tvec_ ) ) );
            sres_   %= trans( expand<E>( eval( tvec_ ) ) );
            osres_  %= trans( expand<E>( eval( tvec_ ) ) );
            refres_ %= trans( expand<E>( eval( trefvec_ ) ) );
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
/*!\brief Testing the conjugate transpose dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the conjugate transpose vector expansion with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the expansion or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testCTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION > 1 )
   {
      using blaze::expand;


      //=====================================================================================
      // Conjugate transpose expansion operation
      //=====================================================================================

      // Conjugate transpose expansion operation with the given vector (runtime)
      {
         test_  = "Conjugate transpose expansion operation with the given vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( expand( vec_, E ) );
            todres_  = ctrans( expand( vec_, E ) );
            tsres_   = ctrans( expand( vec_, E ) );
            tosres_  = ctrans( expand( vec_, E ) );
            trefres_ = ctrans( expand( refvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = ctrans( expand( tvec_, E ) );
            odres_  = ctrans( expand( tvec_, E ) );
            sres_   = ctrans( expand( tvec_, E ) );
            osres_  = ctrans( expand( tvec_, E ) );
            refres_ = ctrans( expand( trefvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose expansion operation with the given vector (compile time)
      {
         test_  = "Conjugate transpose expansion operation with the given vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( expand<E>( vec_ ) );
            todres_  = ctrans( expand<E>( vec_ ) );
            tsres_   = ctrans( expand<E>( vec_ ) );
            tosres_  = ctrans( expand<E>( vec_ ) );
            trefres_ = ctrans( expand<E>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = ctrans( expand<E>( tvec_ ) );
            odres_  = ctrans( expand<E>( tvec_ ) );
            sres_   = ctrans( expand<E>( tvec_ ) );
            osres_  = ctrans( expand<E>( tvec_ ) );
            refres_ = ctrans( expand<E>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose expansion operation with evaluated vector (runtime)
      {
         test_  = "Conjugate transpose expansion operation with evaluated vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( expand( eval( vec_ ), E ) );
            todres_  = ctrans( expand( eval( vec_ ), E ) );
            tsres_   = ctrans( expand( eval( vec_ ), E ) );
            tosres_  = ctrans( expand( eval( vec_ ), E ) );
            trefres_ = ctrans( expand( eval( refvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = ctrans( expand( eval( tvec_ ), E ) );
            odres_  = ctrans( expand( eval( tvec_ ), E ) );
            sres_   = ctrans( expand( eval( tvec_ ), E ) );
            osres_  = ctrans( expand( eval( tvec_ ), E ) );
            refres_ = ctrans( expand( eval( trefvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose expansion operation with evaluated vector (compile time)
      {
         test_  = "Conjugate transpose expansion operation with evaluated vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( expand<E>( eval( vec_ ) ) );
            todres_  = ctrans( expand<E>( eval( vec_ ) ) );
            tsres_   = ctrans( expand<E>( eval( vec_ ) ) );
            tosres_  = ctrans( expand<E>( eval( vec_ ) ) );
            trefres_ = ctrans( expand<E>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   = ctrans( expand<E>( eval( tvec_ ) ) );
            odres_  = ctrans( expand<E>( eval( tvec_ ) ) );
            sres_   = ctrans( expand<E>( eval( tvec_ ) ) );
            osres_  = ctrans( expand<E>( eval( tvec_ ) ) );
            refres_ = ctrans( expand<E>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }


      //=====================================================================================
      // Conjugate transpose expansion with addition assignment
      //=====================================================================================

      // Conjugate transpose expansion with addition assignment with the given vector (runtime)
      {
         test_  = "Conjugate transpose expansion with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += ctrans( expand( vec_, E ) );
            todres_  += ctrans( expand( vec_, E ) );
            tsres_   += ctrans( expand( vec_, E ) );
            tosres_  += ctrans( expand( vec_, E ) );
            trefres_ += ctrans( expand( refvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += ctrans( expand( tvec_, E ) );
            odres_  += ctrans( expand( tvec_, E ) );
            sres_   += ctrans( expand( tvec_, E ) );
            osres_  += ctrans( expand( tvec_, E ) );
            refres_ += ctrans( expand( trefvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose expansion with addition assignment with the given vector (compile time)
      {
         test_  = "Conjugate transpose expansion with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += ctrans( expand<E>( vec_ ) );
            todres_  += ctrans( expand<E>( vec_ ) );
            tsres_   += ctrans( expand<E>( vec_ ) );
            tosres_  += ctrans( expand<E>( vec_ ) );
            trefres_ += ctrans( expand<E>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += ctrans( expand<E>( tvec_ ) );
            odres_  += ctrans( expand<E>( tvec_ ) );
            sres_   += ctrans( expand<E>( tvec_ ) );
            osres_  += ctrans( expand<E>( tvec_ ) );
            refres_ += ctrans( expand<E>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose expansion with addition assignment with evaluated vector (runtime)
      {
         test_  = "Conjugate transpose expansion with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += ctrans( expand( eval( vec_ ), E ) );
            todres_  += ctrans( expand( eval( vec_ ), E ) );
            tsres_   += ctrans( expand( eval( vec_ ), E ) );
            tosres_  += ctrans( expand( eval( vec_ ), E ) );
            trefres_ += ctrans( expand( eval( refvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += ctrans( expand( eval( tvec_ ), E ) );
            odres_  += ctrans( expand( eval( tvec_ ), E ) );
            sres_   += ctrans( expand( eval( tvec_ ), E ) );
            osres_  += ctrans( expand( eval( tvec_ ), E ) );
            refres_ += ctrans( expand( eval( trefvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose expansion with addition assignment with evaluated vector (compile time)
      {
         test_  = "Conjugate transpose expansion with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initTransposeResults();
            tdres_   += ctrans( expand<E>( eval( vec_ ) ) );
            todres_  += ctrans( expand<E>( eval( vec_ ) ) );
            tsres_   += ctrans( expand<E>( eval( vec_ ) ) );
            tosres_  += ctrans( expand<E>( eval( vec_ ) ) );
            trefres_ += ctrans( expand<E>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   += ctrans( expand<E>( eval( tvec_ ) ) );
            odres_  += ctrans( expand<E>( eval( tvec_ ) ) );
            sres_   += ctrans( expand<E>( eval( tvec_ ) ) );
            osres_  += ctrans( expand<E>( eval( tvec_ ) ) );
            refres_ += ctrans( expand<E>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }


      //=====================================================================================
      // Conjugate transpose expansion with subtraction assignment
      //=====================================================================================

      // Conjugate transpose expansion with subtraction assignment with the given vector (runtime)
      {
         test_  = "Conjugate transpose expansion with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= ctrans( expand( vec_, E ) );
            todres_  -= ctrans( expand( vec_, E ) );
            tsres_   -= ctrans( expand( vec_, E ) );
            tosres_  -= ctrans( expand( vec_, E ) );
            trefres_ -= ctrans( expand( refvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= ctrans( expand( tvec_, E ) );
            odres_  -= ctrans( expand( tvec_, E ) );
            sres_   -= ctrans( expand( tvec_, E ) );
            osres_  -= ctrans( expand( tvec_, E ) );
            refres_ -= ctrans( expand( trefvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose expansion with subtraction assignment with the given vector (compile time)
      {
         test_  = "Conjugate transpose expansion with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= ctrans( expand<E>( vec_ ) );
            todres_  -= ctrans( expand<E>( vec_ ) );
            tsres_   -= ctrans( expand<E>( vec_ ) );
            tosres_  -= ctrans( expand<E>( vec_ ) );
            trefres_ -= ctrans( expand<E>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= ctrans( expand<E>( tvec_ ) );
            odres_  -= ctrans( expand<E>( tvec_ ) );
            sres_   -= ctrans( expand<E>( tvec_ ) );
            osres_  -= ctrans( expand<E>( tvec_ ) );
            refres_ -= ctrans( expand<E>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose expansion with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Conjugate transpose expansion with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= ctrans( expand( eval( vec_ ), E ) );
            todres_  -= ctrans( expand( eval( vec_ ), E ) );
            tsres_   -= ctrans( expand( eval( vec_ ), E ) );
            tosres_  -= ctrans( expand( eval( vec_ ), E ) );
            trefres_ -= ctrans( expand( eval( refvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= ctrans( expand( eval( tvec_ ), E ) );
            odres_  -= ctrans( expand( eval( tvec_ ), E ) );
            sres_   -= ctrans( expand( eval( tvec_ ), E ) );
            osres_  -= ctrans( expand( eval( tvec_ ), E ) );
            refres_ -= ctrans( expand( eval( trefvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose expansion with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Conjugate transpose expansion with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initTransposeResults();
            tdres_   -= ctrans( expand<E>( eval( vec_ ) ) );
            todres_  -= ctrans( expand<E>( eval( vec_ ) ) );
            tsres_   -= ctrans( expand<E>( eval( vec_ ) ) );
            tosres_  -= ctrans( expand<E>( eval( vec_ ) ) );
            trefres_ -= ctrans( expand<E>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   -= ctrans( expand<E>( eval( tvec_ ) ) );
            odres_  -= ctrans( expand<E>( eval( tvec_ ) ) );
            sres_   -= ctrans( expand<E>( eval( tvec_ ) ) );
            osres_  -= ctrans( expand<E>( eval( tvec_ ) ) );
            refres_ -= ctrans( expand<E>( eval( trefvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }


      //=====================================================================================
      // Conjugate transpose expansion with Schur product assignment
      //=====================================================================================

      // Conjugate transpose expansion with Schur product assignment with the given vector (runtime)
      {
         test_  = "Conjugate transpose expansion with Schur product assignment with the given vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initTransposeResults();
            tdres_   %= ctrans( expand( vec_, E ) );
            todres_  %= ctrans( expand( vec_, E ) );
            tsres_   %= ctrans( expand( vec_, E ) );
            tosres_  %= ctrans( expand( vec_, E ) );
            trefres_ %= ctrans( expand( refvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   %= ctrans( expand( tvec_, E ) );
            odres_  %= ctrans( expand( tvec_, E ) );
            sres_   %= ctrans( expand( tvec_, E ) );
            osres_  %= ctrans( expand( tvec_, E ) );
            refres_ %= ctrans( expand( trefvec_, E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose expansion with Schur product assignment with the given vector (compile time)
      {
         test_  = "Conjugate transpose expansion with Schur product assignment with the given vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initTransposeResults();
            tdres_   %= ctrans( expand<E>( vec_ ) );
            todres_  %= ctrans( expand<E>( vec_ ) );
            tsres_   %= ctrans( expand<E>( vec_ ) );
            tosres_  %= ctrans( expand<E>( vec_ ) );
            trefres_ %= ctrans( expand<E>( refvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   %= ctrans( expand<E>( tvec_ ) );
            odres_  %= ctrans( expand<E>( tvec_ ) );
            sres_   %= ctrans( expand<E>( tvec_ ) );
            osres_  %= ctrans( expand<E>( tvec_ ) );
            refres_ %= ctrans( expand<E>( trefvec_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose expansion with Schur product assignment with evaluated vector (runtime)
      {
         test_  = "Conjugate transpose expansion with Schur product assignment with evaluated vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initTransposeResults();
            tdres_   %= ctrans( expand( eval( vec_ ), E ) );
            todres_  %= ctrans( expand( eval( vec_ ), E ) );
            tsres_   %= ctrans( expand( eval( vec_ ), E ) );
            tosres_  %= ctrans( expand( eval( vec_ ), E ) );
            trefres_ %= ctrans( expand( eval( refvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   %= ctrans( expand( eval( tvec_ ), E ) );
            odres_  %= ctrans( expand( eval( tvec_ ), E ) );
            sres_   %= ctrans( expand( eval( tvec_ ), E ) );
            osres_  %= ctrans( expand( eval( tvec_ ), E ) );
            refres_ %= ctrans( expand( eval( trefvec_ ), E ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkResults<TVT>();
      }

      // Conjugate transpose expansion with Schur product assignment with evaluated vector (compile time)
      {
         test_  = "Conjugate transpose expansion with Schur product assignment with evaluated vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initTransposeResults();
            tdres_   %= ctrans( expand<E>( eval( vec_ ) ) );
            todres_  %= ctrans( expand<E>( eval( vec_ ) ) );
            tsres_   %= ctrans( expand<E>( eval( vec_ ) ) );
            tosres_  %= ctrans( expand<E>( eval( vec_ ) ) );
            trefres_ %= ctrans( expand<E>( eval( refvec_ ) ) );
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkTransposeResults<VT>();

         try {
            initResults();
            dres_   %= ctrans( expand<E>( eval( tvec_ ) ) );
            odres_  %= ctrans( expand<E>( eval( tvec_ ) ) );
            sres_   %= ctrans( expand<E>( eval( tvec_ ) ) );
            osres_  %= ctrans( expand<E>( eval( tvec_ ) ) );
            refres_ %= ctrans( expand<E>( eval( trefvec_ ) ) );
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
/*!\brief Testing the abs dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the abs vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// expansion or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testAbsOperation()
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
/*!\brief Testing the conjugate dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the conjugate vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// expansion or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testConjOperation()
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
/*!\brief Testing the \a real dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the \a real vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// expansion or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testRealOperation()
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
/*!\brief Testing the \a imag dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the \a imag vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// expansion or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testImagOperation()
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
/*!\brief Testing the evaluated dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the evalated vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// expansion or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testEvalOperation()
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
/*!\brief Testing the serialized dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the serialized vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// expansion or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testSerialOperation()
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
/*!\brief Testing the non-aliased dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the non-aliased vector expansion with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the expansion or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testNoAliasOperation()
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
/*!\brief Testing the non-SIMD dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the non-SIMD vector expansion with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the expansion or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testNoSIMDOperation()
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
/*!\brief Testing the submatrix-wise dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the submatrix-wise vector expansion with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the addition or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testSubmatrixOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION > 1 )
   {
      using blaze::expand;

      if( vec_.size() == 0UL || E == 0UL )
         return;


      //=====================================================================================
      // Submatrix-wise expansion
      //=====================================================================================

      // Submatrix-wise expansion with the given vector (runtime)
      {
         test_  = "Submatrix-wise expansion with the given vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( expand( refvec_, E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) = submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) = submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) = submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) = submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) = submatrix( expand( trefvec_, E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Submatrix-wise expansion with the given vector (compile time)
      {
         test_  = "Submatrix-wise expansion with the given vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( expand<E>( refvec_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) = submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) = submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) = submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) = submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) = submatrix( expand<E>( trefvec_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Submatrix-wise expansion with evaluated vector (runtime)
      {
         test_  = "Submatrix-wise expansion with evaluated vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( expand( eval( refvec_ ), E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) = submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) = submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) = submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) = submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) = submatrix( expand( eval( trefvec_ ), E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Submatrix-wise expansion with evaluated vector (compile time)
      {
         test_  = "Submatrix-wise expansion with evaluated vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) = submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) = submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) = submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) = submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) = submatrix( expand<E>( eval( refvec_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) = submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) = submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) = submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) = submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) = submatrix( expand<E>( eval( trefvec_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Submatrix-wise expansion with addition assignment
      //=====================================================================================

      // Submatrix-wise expansion with addition assignment with the given vector (runtime)
      {
         test_  = "Submatrix-wise expansion with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( expand( refvec_, E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) += submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) += submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) += submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) += submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) += submatrix( expand( trefvec_, E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Submatrix-wise expansion with addition assignment with the given vector (compile time)
      {
         test_  = "Submatrix-wise expansion with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( expand<E>( refvec_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) += submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) += submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) += submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) += submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) += submatrix( expand<E>( trefvec_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Submatrix-wise expansion with addition assignment with evaluated vector (runtime)
      {
         test_  = "Submatrix-wise expansion with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( expand( eval( refvec_ ), E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) += submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) += submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) += submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) += submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) += submatrix( expand( eval( trefvec_ ), E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Submatrix-wise expansion with addition assignment with evaluated vector (compile time)
      {
         test_  = "Submatrix-wise expansion with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) += submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) += submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) += submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) += submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) += submatrix( expand<E>( eval( refvec_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) += submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) += submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) += submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) += submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) += submatrix( expand<E>( eval( trefvec_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Submatrix-wise expansion with subtraction assignment
      //=====================================================================================

      // Submatrix-wise expansion with subtraction assignment with the given vector (runtime)
      {
         test_  = "Submatrix-wise expansion with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( expand( refvec_, E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) -= submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) -= submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) -= submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) -= submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) -= submatrix( expand( trefvec_, E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Submatrix-wise expansion with subtraction assignment with the given vector (compile time)
      {
         test_  = "Submatrix-wise expansion with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( expand<E>( refvec_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) -= submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) -= submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) -= submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) -= submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) -= submatrix( expand<E>( trefvec_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Submatrix-wise expansion with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Submatrix-wise expansion with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( expand( eval( refvec_ ), E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) -= submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) -= submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) -= submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) -= submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) -= submatrix( expand( eval( trefvec_ ), E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Submatrix-wise expansion with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Submatrix-wise expansion with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) -= submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) -= submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) -= submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) -= submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) -= submatrix( expand<E>( eval( refvec_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) -= submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) -= submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) -= submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) -= submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) -= submatrix( expand<E>( eval( trefvec_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Submatrix-wise expansion with Schur product assignment
      //=====================================================================================

      // Submatrix-wise expansion with Schur product assignment with the given vector (runtime)
      {
         test_  = "Submatrix-wise expansion with Schur product assignment with the given vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( expand( vec_, E )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( expand( refvec_, E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) %= submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) %= submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) %= submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) %= submatrix( expand( tvec_, E )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) %= submatrix( expand( trefvec_, E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Submatrix-wise expansion with Schur product assignment with the given vector (compile time)
      {
         test_  = "Submatrix-wise expansion with Schur product assignment with the given vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( expand<E>( vec_ )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( expand<E>( refvec_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) %= submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) %= submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) %= submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) %= submatrix( expand<E>( tvec_ )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) %= submatrix( expand<E>( trefvec_ ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Submatrix-wise expansion with Schur product assignment with evaluated vector (runtime)
      {
         test_  = "Submatrix-wise expansion with Schur product assignment with evaluated vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( expand( eval( vec_ ), E )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( expand( eval( refvec_ ), E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) %= submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) %= submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) %= submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) %= submatrix( expand( eval( tvec_ ), E )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) %= submatrix( expand( eval( trefvec_ ), E ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Submatrix-wise expansion with Schur product assignment with evaluated vector (compile time)
      {
         test_  = "Submatrix-wise expansion with Schur product assignment with evaluated vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t row=0UL, m=0UL; row<vec_.size(); row+=m ) {
               m = blaze::rand<size_t>( 1UL, vec_.size() - row );
               for( size_t column=0UL, n=0UL; column<E; column+=n ) {
                  n = blaze::rand<size_t>( 1UL, E - column );
                  submatrix( dres_  , row, column, m, n ) %= submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( odres_ , row, column, m, n ) %= submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( sres_  , row, column, m, n ) %= submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( osres_ , row, column, m, n ) %= submatrix( expand<E>( eval( vec_ ) )   , row, column, m, n );
                  submatrix( refres_, row, column, m, n ) %= submatrix( expand<E>( eval( refvec_ ) ), row, column, m, n );
               }
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t row=0UL, m=0UL; row<E; row+=m ) {
               m = blaze::rand<size_t>( 1UL, E - row );
               for( size_t column=0UL, n=0UL; column<tvec_.size(); column+=n ) {
                  n = blaze::rand<size_t>( 1UL, tvec_.size() - column );
                  submatrix( tdres_  , row, column, m, n ) %= submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( todres_ , row, column, m, n ) %= submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( tsres_  , row, column, m, n ) %= submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( tosres_ , row, column, m, n ) %= submatrix( expand<E>( eval( tvec_ ) )   , row, column, m, n );
                  submatrix( trefres_, row, column, m, n ) %= submatrix( expand<E>( eval( trefvec_ ) ), row, column, m, n );
               }
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

      // Out-of-bounds submatrix construction (invalid number of rows)
      {
         test_  = "Out-of-bounds submatrix construction (invalid number of rows)";
         error_ = "Setup of out-of-bounds submatrix succeeded";

         try {
            auto sm = submatrix( expand( vec_, E ), 1UL, 0UL, vec_.size(), E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( expand<E>( vec_ ), 1UL, 0UL, vec_.size(), E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( expand( tvec_, E ), 1UL, 0UL, E, tvec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( expand<E>( tvec_ ), 1UL, 0UL, E, tvec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
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
            auto sm = submatrix( expand( vec_, E ), 0UL, 1UL, vec_.size(), E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( expand<E>( vec_ ), 0UL, 1UL, vec_.size(), E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( expand( tvec_, E ), 0UL, 1UL, E, tvec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( expand<E>( tvec_ ), 0UL, 1UL, E, tvec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
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
            auto sm = submatrix( expand( vec_, E ), vec_.size(), 0UL, 1UL, E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( expand<E>( vec_ ), vec_.size(), 0UL, 1UL, E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( expand( tvec_, E ), E, 0UL, 1UL, tvec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( expand<E>( tvec_ ), E, 0UL, 1UL, tvec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
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
            auto sm = submatrix( expand( vec_, E ), 0UL, E, vec_.size(), 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( expand<E>( vec_ ), 0UL, E, vec_.size(), 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( expand( tvec_, E ), 0UL, tvec_.size(), E, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid submatrix specification" );
         }

         try {
            auto sm = submatrix( expand<E>( tvec_ ), 0UL, tvec_.size(), E, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
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
/*!\brief Skipping the submatrix-wise dense vector expansion operation.
//
// \return void
//
// This function is called in case the submatrix-wise dense vector expansion operation is not
// available for the given vector type \a VT.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testSubmatrixOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the row-wise dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the row-wise vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// addition or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testRowOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_ROW_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ROW_OPERATION > 1 )
   {
      using blaze::expand;

      if( vec_.size() == 0UL || E == 0UL )
         return;


      //=====================================================================================
      // Row-wise expansion
      //=====================================================================================

      // Row-wise expansion with the given vector (runtime)
      {
         test_  = "Row-wise expansion with the given vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) = row( expand( vec_, E ), i );
               row( odres_ , i ) = row( expand( vec_, E ), i );
               row( sres_  , i ) = row( expand( vec_, E ), i );
               row( osres_ , i ) = row( expand( vec_, E ), i );
               row( refres_, i ) = row( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) = row( expand( tvec_, E ), i );
               row( todres_ , i ) = row( expand( tvec_, E ), i );
               row( tsres_  , i ) = row( expand( tvec_, E ), i );
               row( tosres_ , i ) = row( expand( tvec_, E ), i );
               row( trefres_, i ) = row( expand( trefvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Row-wise expansion with the given vector (compile time)
      {
         test_  = "Row-wise expansion with the given vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) = row( expand<E>( vec_ ), i );
               row( odres_ , i ) = row( expand<E>( vec_ ), i );
               row( sres_  , i ) = row( expand<E>( vec_ ), i );
               row( osres_ , i ) = row( expand<E>( vec_ ), i );
               row( refres_, i ) = row( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) = row( expand<E>( tvec_ ), i );
               row( todres_ , i ) = row( expand<E>( tvec_ ), i );
               row( tsres_  , i ) = row( expand<E>( tvec_ ), i );
               row( tosres_ , i ) = row( expand<E>( tvec_ ), i );
               row( trefres_, i ) = row( expand<E>( trefvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Row-wise expansion with evaluated vector (runtime)
      {
         test_  = "Row-wise expansion with evaluated vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) = row( expand( eval( vec_ ), E ), i );
               row( odres_ , i ) = row( expand( eval( vec_ ), E ), i );
               row( sres_  , i ) = row( expand( eval( vec_ ), E ), i );
               row( osres_ , i ) = row( expand( eval( vec_ ), E ), i );
               row( refres_, i ) = row( expand( eval( refvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) = row( expand( eval( tvec_ ), E ), i );
               row( todres_ , i ) = row( expand( eval( tvec_ ), E ), i );
               row( tsres_  , i ) = row( expand( eval( tvec_ ), E ), i );
               row( tosres_ , i ) = row( expand( eval( tvec_ ), E ), i );
               row( trefres_, i ) = row( expand( eval( trefvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Row-wise expansion with evaluated vector (compile time)
      {
         test_  = "Row-wise expansion with evaluated vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) = row( expand<E>( eval( vec_ ) ), i );
               row( odres_ , i ) = row( expand<E>( eval( vec_ ) ), i );
               row( sres_  , i ) = row( expand<E>( eval( vec_ ) ), i );
               row( osres_ , i ) = row( expand<E>( eval( vec_ ) ), i );
               row( refres_, i ) = row( expand<E>( eval( refvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) = row( expand<E>( eval( tvec_ ) ), i );
               row( todres_ , i ) = row( expand<E>( eval( tvec_ ) ), i );
               row( tsres_  , i ) = row( expand<E>( eval( tvec_ ) ), i );
               row( tosres_ , i ) = row( expand<E>( eval( tvec_ ) ), i );
               row( trefres_, i ) = row( expand<E>( eval( trefvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Row-wise expansion with addition assignment
      //=====================================================================================

      // Row-wise expansion with addition assignment with the given vector (runtime)
      {
         test_  = "Row-wise expansion with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) += row( expand( vec_, E ), i );
               row( odres_ , i ) += row( expand( vec_, E ), i );
               row( sres_  , i ) += row( expand( vec_, E ), i );
               row( osres_ , i ) += row( expand( vec_, E ), i );
               row( refres_, i ) += row( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) += row( expand( tvec_, E ), i );
               row( todres_ , i ) += row( expand( tvec_, E ), i );
               row( tsres_  , i ) += row( expand( tvec_, E ), i );
               row( tosres_ , i ) += row( expand( tvec_, E ), i );
               row( trefres_, i ) += row( expand( trefvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Row-wise expansion with addition assignment with the given vector (compile time)
      {
         test_  = "Row-wise expansion with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) += row( expand<E>( vec_ ), i );
               row( odres_ , i ) += row( expand<E>( vec_ ), i );
               row( sres_  , i ) += row( expand<E>( vec_ ), i );
               row( osres_ , i ) += row( expand<E>( vec_ ), i );
               row( refres_, i ) += row( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) += row( expand<E>( tvec_ ), i );
               row( todres_ , i ) += row( expand<E>( tvec_ ), i );
               row( tsres_  , i ) += row( expand<E>( tvec_ ), i );
               row( tosres_ , i ) += row( expand<E>( tvec_ ), i );
               row( trefres_, i ) += row( expand<E>( trefvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Row-wise expansion with addition assignment with evaluated vector (runtime)
      {
         test_  = "Row-wise expansion with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) += row( expand( eval( vec_ ), E ), i );
               row( odres_ , i ) += row( expand( eval( vec_ ), E ), i );
               row( sres_  , i ) += row( expand( eval( vec_ ), E ), i );
               row( osres_ , i ) += row( expand( eval( vec_ ), E ), i );
               row( refres_, i ) += row( expand( eval( refvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) += row( expand( eval( tvec_ ), E ), i );
               row( todres_ , i ) += row( expand( eval( tvec_ ), E ), i );
               row( tsres_  , i ) += row( expand( eval( tvec_ ), E ), i );
               row( tosres_ , i ) += row( expand( eval( tvec_ ), E ), i );
               row( trefres_, i ) += row( expand( eval( trefvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Row-wise expansion with addition assignment with evaluated vector (compile time)
      {
         test_  = "Row-wise expansion with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) += row( expand<E>( eval( vec_ ) ), i );
               row( odres_ , i ) += row( expand<E>( eval( vec_ ) ), i );
               row( sres_  , i ) += row( expand<E>( eval( vec_ ) ), i );
               row( osres_ , i ) += row( expand<E>( eval( vec_ ) ), i );
               row( refres_, i ) += row( expand<E>( eval( refvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) += row( expand<E>( eval( tvec_ ) ), i );
               row( todres_ , i ) += row( expand<E>( eval( tvec_ ) ), i );
               row( tsres_  , i ) += row( expand<E>( eval( tvec_ ) ), i );
               row( tosres_ , i ) += row( expand<E>( eval( tvec_ ) ), i );
               row( trefres_, i ) += row( expand<E>( eval( trefvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Row-wise expansion with subtraction assignment
      //=====================================================================================

      // Row-wise expansion with subtraction assignment with the given vector (runtime)
      {
         test_  = "Row-wise expansion with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) -= row( expand( vec_, E ), i );
               row( odres_ , i ) -= row( expand( vec_, E ), i );
               row( sres_  , i ) -= row( expand( vec_, E ), i );
               row( osres_ , i ) -= row( expand( vec_, E ), i );
               row( refres_, i ) -= row( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) -= row( expand( tvec_, E ), i );
               row( todres_ , i ) -= row( expand( tvec_, E ), i );
               row( tsres_  , i ) -= row( expand( tvec_, E ), i );
               row( tosres_ , i ) -= row( expand( tvec_, E ), i );
               row( trefres_, i ) -= row( expand( trefvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Row-wise expansion with subtraction assignment with the given vector (compile time)
      {
         test_  = "Row-wise expansion with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) -= row( expand<E>( vec_ ), i );
               row( odres_ , i ) -= row( expand<E>( vec_ ), i );
               row( sres_  , i ) -= row( expand<E>( vec_ ), i );
               row( osres_ , i ) -= row( expand<E>( vec_ ), i );
               row( refres_, i ) -= row( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) -= row( expand<E>( tvec_ ), i );
               row( todres_ , i ) -= row( expand<E>( tvec_ ), i );
               row( tsres_  , i ) -= row( expand<E>( tvec_ ), i );
               row( tosres_ , i ) -= row( expand<E>( tvec_ ), i );
               row( trefres_, i ) -= row( expand<E>( trefvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Row-wise expansion with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Row-wise expansion with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) -= row( expand( eval( vec_ ), E ), i );
               row( odres_ , i ) -= row( expand( eval( vec_ ), E ), i );
               row( sres_  , i ) -= row( expand( eval( vec_ ), E ), i );
               row( osres_ , i ) -= row( expand( eval( vec_ ), E ), i );
               row( refres_, i ) -= row( expand( eval( refvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) -= row( expand( eval( tvec_ ), E ), i );
               row( todres_ , i ) -= row( expand( eval( tvec_ ), E ), i );
               row( tsres_  , i ) -= row( expand( eval( tvec_ ), E ), i );
               row( tosres_ , i ) -= row( expand( eval( tvec_ ), E ), i );
               row( trefres_, i ) -= row( expand( eval( trefvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Row-wise expansion with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Row-wise expansion with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) -= row( expand<E>( eval( vec_ ) ), i );
               row( odres_ , i ) -= row( expand<E>( eval( vec_ ) ), i );
               row( sres_  , i ) -= row( expand<E>( eval( vec_ ) ), i );
               row( osres_ , i ) -= row( expand<E>( eval( vec_ ) ), i );
               row( refres_, i ) -= row( expand<E>( eval( refvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) -= row( expand<E>( eval( tvec_ ) ), i );
               row( todres_ , i ) -= row( expand<E>( eval( tvec_ ) ), i );
               row( tsres_  , i ) -= row( expand<E>( eval( tvec_ ) ), i );
               row( tosres_ , i ) -= row( expand<E>( eval( tvec_ ) ), i );
               row( trefres_, i ) -= row( expand<E>( eval( trefvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Row-wise expansion with multiplication assignment
      //=====================================================================================

      // Row-wise expansion with multiplication assignment with the given vector (runtime)
      {
         test_  = "Row-wise expansion with multiplication assignment with the given vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) *= row( expand( vec_, E ), i );
               row( odres_ , i ) *= row( expand( vec_, E ), i );
               row( sres_  , i ) *= row( expand( vec_, E ), i );
               row( osres_ , i ) *= row( expand( vec_, E ), i );
               row( refres_, i ) *= row( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) *= row( expand( tvec_, E ), i );
               row( todres_ , i ) *= row( expand( tvec_, E ), i );
               row( tsres_  , i ) *= row( expand( tvec_, E ), i );
               row( tosres_ , i ) *= row( expand( tvec_, E ), i );
               row( trefres_, i ) *= row( expand( trefvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Row-wise expansion with multiplication assignment with the given vector (compile time)
      {
         test_  = "Row-wise expansion with multiplication assignment with the given vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) *= row( expand<E>( vec_ ), i );
               row( odres_ , i ) *= row( expand<E>( vec_ ), i );
               row( sres_  , i ) *= row( expand<E>( vec_ ), i );
               row( osres_ , i ) *= row( expand<E>( vec_ ), i );
               row( refres_, i ) *= row( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) *= row( expand<E>( tvec_ ), i );
               row( todres_ , i ) *= row( expand<E>( tvec_ ), i );
               row( tsres_  , i ) *= row( expand<E>( tvec_ ), i );
               row( tosres_ , i ) *= row( expand<E>( tvec_ ), i );
               row( trefres_, i ) *= row( expand<E>( trefvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Row-wise expansion with multiplication assignment with evaluated vector (runtime)
      {
         test_  = "Row-wise expansion with multiplication assignment with evaluated vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) *= row( expand( eval( vec_ ), E ), i );
               row( odres_ , i ) *= row( expand( eval( vec_ ), E ), i );
               row( sres_  , i ) *= row( expand( eval( vec_ ), E ), i );
               row( osres_ , i ) *= row( expand( eval( vec_ ), E ), i );
               row( refres_, i ) *= row( expand( eval( refvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) *= row( expand( eval( tvec_ ), E ), i );
               row( todres_ , i ) *= row( expand( eval( tvec_ ), E ), i );
               row( tsres_  , i ) *= row( expand( eval( tvec_ ), E ), i );
               row( tosres_ , i ) *= row( expand( eval( tvec_ ), E ), i );
               row( trefres_, i ) *= row( expand( eval( trefvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Row-wise expansion with multiplication assignment with evaluated vector (compile time)
      {
         test_  = "Row-wise expansion with multiplication assignment with evaluated vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               row( dres_  , i ) *= row( expand<E>( eval( vec_ ) ), i );
               row( odres_ , i ) *= row( expand<E>( eval( vec_ ) ), i );
               row( sres_  , i ) *= row( expand<E>( eval( vec_ ) ), i );
               row( osres_ , i ) *= row( expand<E>( eval( vec_ ) ), i );
               row( refres_, i ) *= row( expand<E>( eval( refvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<E; ++i ) {
               row( tdres_  , i ) *= row( expand<E>( eval( tvec_ ) ), i );
               row( todres_ , i ) *= row( expand<E>( eval( tvec_ ) ), i );
               row( tsres_  , i ) *= row( expand<E>( eval( tvec_ ) ), i );
               row( tosres_ , i ) *= row( expand<E>( eval( tvec_ ) ), i );
               row( trefres_, i ) *= row( expand<E>( eval( trefvec_ ) ), i );
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

      // Out-of-bounds access (invalid row index)
      {
         test_  = "Out-of-bounds row construction (invalid row index)";
         error_ = "Setup of out-of-bounds row succeeded";

         try {
            auto r = row( expand( vec_, E ), vec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = row( expand<E>( vec_ ), vec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = row( expand( tvec_, E ), E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = row( expand<E>( tvec_ ), E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
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
/*!\brief Skipping the row-wise dense vector expansion operation.
//
// \return void
//
// This function is called in case the row-wise dense vector expansion operation is not
// available for the given vector type \a VT.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testRowOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the rows-wise dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the rows-wise vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// addition or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testRowsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_ROWS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ROWS_OPERATION > 1 )
   {
      using blaze::expand;

      if( vec_.size() == 0UL || E == 0UL )
         return;


      std::vector<size_t> indices( vec_.size() );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );

      std::vector<size_t> tindices( E );
      std::iota( tindices.begin(), tindices.end(), 0UL );
      std::random_shuffle( tindices.begin(), tindices.end() );


      //=====================================================================================
      // Rows-wise expansion
      //=====================================================================================

      // Rows-wise expansion with the given vector (runtime)
      {
         test_  = "Rows-wise expansion with the given vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( expand( vec_, E ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( expand( vec_, E ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( expand( vec_, E ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( expand( vec_, E ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( expand( refvec_, E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) = rows( expand( tvec_, E ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) = rows( expand( tvec_, E ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) = rows( expand( tvec_, E ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) = rows( expand( tvec_, E ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) = rows( expand( trefvec_, E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Rows-wise expansion with the given vector (compile time)
      {
         test_  = "Rows-wise expansion with the given vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( expand<E>( vec_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( expand<E>( vec_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( expand<E>( vec_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( expand<E>( vec_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( expand<E>( refvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) = rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) = rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) = rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) = rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) = rows( expand<E>( trefvec_ ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Rows-wise expansion with evaluated vector (runtime)
      {
         test_  = "Rows-wise expansion with evaluated vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( expand( eval( refvec_ ), E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) = rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) = rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) = rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) = rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) = rows( expand( eval( trefvec_ ), E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Rows-wise expansion with evaluated vector (compile time)
      {
         test_  = "Rows-wise expansion with evaluated vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) = rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) = rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) = rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) = rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) = rows( expand<E>( eval( refvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) = rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) = rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) = rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) = rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) = rows( expand<E>( eval( trefvec_ ) ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Rows-wise expansion with addition assignment
      //=====================================================================================

      // Rows-wise expansion with addition assignment with the given vector (runtime)
      {
         test_  = "Rows-wise expansion with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( expand( vec_, E ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( expand( vec_, E ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( expand( vec_, E ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( expand( vec_, E ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( expand( refvec_, E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) += rows( expand( tvec_, E ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) += rows( expand( tvec_, E ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) += rows( expand( tvec_, E ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) += rows( expand( tvec_, E ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) += rows( expand( trefvec_, E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Rows-wise expansion with addition assignment with the given vector (compile time)
      {
         test_  = "Rows-wise expansion with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( expand<E>( vec_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( expand<E>( vec_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( expand<E>( vec_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( expand<E>( vec_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( expand<E>( refvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) += rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) += rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) += rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) += rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) += rows( expand<E>( trefvec_ ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Rows-wise expansion with addition assignment with evaluated vector (runtime)
      {
         test_  = "Rows-wise expansion with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( expand( eval( refvec_ ), E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) += rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) += rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) += rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) += rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) += rows( expand( eval( trefvec_ ), E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Rows-wise expansion with addition assignment with evaluated vector (compile time)
      {
         test_  = "Rows-wise expansion with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) += rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) += rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) += rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) += rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) += rows( expand<E>( eval( refvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) += rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) += rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) += rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) += rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) += rows( expand<E>( eval( trefvec_ ) ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Rows-wise expansion with subtraction assignment
      //=====================================================================================

      // Rows-wise expansion with subtraction assignment with the given vector (runtime)
      {
         test_  = "Rows-wise expansion with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( expand( vec_, E ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( expand( vec_, E ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( expand( vec_, E ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( expand( vec_, E ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( expand( refvec_, E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) -= rows( expand( tvec_, E ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) -= rows( expand( tvec_, E ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) -= rows( expand( tvec_, E ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) -= rows( expand( tvec_, E ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) -= rows( expand( trefvec_, E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Rows-wise expansion with subtraction assignment with the given vector (compile time)
      {
         test_  = "Rows-wise expansion with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( expand<E>( vec_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( expand<E>( vec_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( expand<E>( vec_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( expand<E>( vec_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( expand<E>( refvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) -= rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) -= rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) -= rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) -= rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) -= rows( expand<E>( trefvec_ ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Rows-wise expansion with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Rows-wise expansion with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( expand( eval( refvec_ ), E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) -= rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) -= rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) -= rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) -= rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) -= rows( expand( eval( trefvec_ ), E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Rows-wise expansion with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Rows-wise expansion with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) -= rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) -= rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) -= rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) -= rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) -= rows( expand<E>( eval( refvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) -= rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) -= rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) -= rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) -= rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) -= rows( expand<E>( eval( trefvec_ ) ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Rows-wise expansion with Schur product assignment
      //=====================================================================================

      // Rows-wise expansion with Schur product assignment with the given vector (runtime)
      {
         test_  = "Rows-wise expansion with Schur product assignment with the given vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( expand( vec_, E ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( expand( vec_, E ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( expand( vec_, E ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( expand( vec_, E ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( expand( refvec_, E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) %= rows( expand( tvec_, E ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) %= rows( expand( tvec_, E ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) %= rows( expand( tvec_, E ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) %= rows( expand( tvec_, E ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) %= rows( expand( trefvec_, E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Rows-wise expansion with Schur product assignment with the given vector (compile time)
      {
         test_  = "Rows-wise expansion with Schur product assignment with the given vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( expand<E>( vec_ ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( expand<E>( vec_ ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( expand<E>( vec_ ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( expand<E>( vec_ ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( expand<E>( refvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) %= rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) %= rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) %= rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) %= rows( expand<E>( tvec_ ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) %= rows( expand<E>( trefvec_ ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Rows-wise expansion with Schur product assignment with evaluated vector (runtime)
      {
         test_  = "Rows-wise expansion with Schur product assignment with evaluated vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( expand( eval( vec_ ), E ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( expand( eval( refvec_ ), E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) %= rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) %= rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) %= rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) %= rows( expand( eval( tvec_ ), E ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) %= rows( expand( eval( trefvec_ ), E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Rows-wise expansion with Schur product assignment with evaluated vector (compile time)
      {
         test_  = "Rows-wise expansion with Schur product assignment with evaluated vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               rows( dres_  , &indices[index], n ) %= rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( odres_ , &indices[index], n ) %= rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( sres_  , &indices[index], n ) %= rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( osres_ , &indices[index], n ) %= rows( expand<E>( eval( vec_ ) ), &indices[index], n );
               rows( refres_, &indices[index], n ) %= rows( expand<E>( eval( refvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               rows( tdres_  , &tindices[index], n ) %= rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( todres_ , &tindices[index], n ) %= rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( tsres_  , &tindices[index], n ) %= rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( tosres_ , &tindices[index], n ) %= rows( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               rows( trefres_, &tindices[index], n ) %= rows( expand<E>( eval( trefvec_ ) ), &tindices[index], n );
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

      // Out-of-bounds access (invalid row index; initializer_list)
      {
         test_  = "Out-of-bounds row selection construction (invalid row index; initializer_list)";
         error_ = "Setup of out-of-bounds row selection succeeded";

         try {
            auto r = rows( expand( vec_, E ), { vec_.size() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( expand<E>( vec_ ), { vec_.size() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( expand( tvec_, E ), { E } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( expand<E>( tvec_ ), { E } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
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
            auto r = rows( expand( vec_, E ), [index=vec_.size()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( expand<E>( vec_ ), [index=vec_.size()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( expand( tvec_, E ), [index=E]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
                << "   Result:\n" << r << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid row access index" );
         }

         try {
            auto r = rows( expand<E>( tvec_ ), [index=E]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
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
/*!\brief Skipping the rows-wise dense vector expansion operation.
//
// \return void
//
// This function is called in case the rows-wise dense vector expansion operation is not
// available for the given vector type \a VT.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testRowsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the column-wise dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the column-wise vector expansion with plain assignment, addition
// assignment, subtraction assignment, and Schur product assignment. In case any error resulting
// from the addition or the subsequent assignment is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testColumnOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION
   if( BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION > 1 )
   {
      using blaze::expand;

      if( vec_.size() == 0UL || E == 0UL )
         return;


      //=====================================================================================
      // Column-wise expansion
      //=====================================================================================

      // Column-wise expansion with the given vector (runtime)
      {
         test_  = "Column-wise expansion with the given vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) = column( expand( vec_, E ), i );
               column( odres_ , i ) = column( expand( vec_, E ), i );
               column( sres_  , i ) = column( expand( vec_, E ), i );
               column( osres_ , i ) = column( expand( vec_, E ), i );
               column( refres_, i ) = column( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) = column( expand( tvec_, E ), i );
               column( todres_ , i ) = column( expand( tvec_, E ), i );
               column( tsres_  , i ) = column( expand( tvec_, E ), i );
               column( tosres_ , i ) = column( expand( tvec_, E ), i );
               column( trefres_, i ) = column( expand( trefvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Column-wise expansion with the given vector (compile time)
      {
         test_  = "Column-wise expansion with the given vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) = column( expand<E>( vec_ ), i );
               column( odres_ , i ) = column( expand<E>( vec_ ), i );
               column( sres_  , i ) = column( expand<E>( vec_ ), i );
               column( osres_ , i ) = column( expand<E>( vec_ ), i );
               column( refres_, i ) = column( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) = column( expand<E>( tvec_ ), i );
               column( todres_ , i ) = column( expand<E>( tvec_ ), i );
               column( tsres_  , i ) = column( expand<E>( tvec_ ), i );
               column( tosres_ , i ) = column( expand<E>( tvec_ ), i );
               column( trefres_, i ) = column( expand<E>( trefvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Column-wise expansion with evaluated vector (runtime)
      {
         test_  = "Column-wise expansion with evaluated vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) = column( expand( eval( vec_ ), E ), i );
               column( odres_ , i ) = column( expand( eval( vec_ ), E ), i );
               column( sres_  , i ) = column( expand( eval( vec_ ), E ), i );
               column( osres_ , i ) = column( expand( eval( vec_ ), E ), i );
               column( refres_, i ) = column( expand( eval( refvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) = column( expand( eval( tvec_ ), E ), i );
               column( todres_ , i ) = column( expand( eval( tvec_ ), E ), i );
               column( tsres_  , i ) = column( expand( eval( tvec_ ), E ), i );
               column( tosres_ , i ) = column( expand( eval( tvec_ ), E ), i );
               column( trefres_, i ) = column( expand( eval( trefvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Column-wise expansion with evaluated vector (compile time)
      {
         test_  = "Column-wise expansion with evaluated vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) = column( expand<E>( eval( vec_ ) ), i );
               column( odres_ , i ) = column( expand<E>( eval( vec_ ) ), i );
               column( sres_  , i ) = column( expand<E>( eval( vec_ ) ), i );
               column( osres_ , i ) = column( expand<E>( eval( vec_ ) ), i );
               column( refres_, i ) = column( expand<E>( eval( refvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) = column( expand<E>( eval( tvec_ ) ), i );
               column( todres_ , i ) = column( expand<E>( eval( tvec_ ) ), i );
               column( tsres_  , i ) = column( expand<E>( eval( tvec_ ) ), i );
               column( tosres_ , i ) = column( expand<E>( eval( tvec_ ) ), i );
               column( trefres_, i ) = column( expand<E>( eval( trefvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Column-wise expansion with addition assignment
      //=====================================================================================

      // Column-wise expansion with addition assignment with the given vector (runtime)
      {
         test_  = "Column-wise expansion with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) += column( expand( vec_, E ), i );
               column( odres_ , i ) += column( expand( vec_, E ), i );
               column( sres_  , i ) += column( expand( vec_, E ), i );
               column( osres_ , i ) += column( expand( vec_, E ), i );
               column( refres_, i ) += column( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) += column( expand( tvec_, E ), i );
               column( todres_ , i ) += column( expand( tvec_, E ), i );
               column( tsres_  , i ) += column( expand( tvec_, E ), i );
               column( tosres_ , i ) += column( expand( tvec_, E ), i );
               column( trefres_, i ) += column( expand( trefvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Column-wise expansion with addition assignment with the given vector (compile time)
      {
         test_  = "Column-wise expansion with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) += column( expand<E>( vec_ ), i );
               column( odres_ , i ) += column( expand<E>( vec_ ), i );
               column( sres_  , i ) += column( expand<E>( vec_ ), i );
               column( osres_ , i ) += column( expand<E>( vec_ ), i );
               column( refres_, i ) += column( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) += column( expand<E>( tvec_ ), i );
               column( todres_ , i ) += column( expand<E>( tvec_ ), i );
               column( tsres_  , i ) += column( expand<E>( tvec_ ), i );
               column( tosres_ , i ) += column( expand<E>( tvec_ ), i );
               column( trefres_, i ) += column( expand<E>( trefvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Column-wise expansion with addition assignment with evaluated vector (runtime)
      {
         test_  = "Column-wise expansion with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) += column( expand( eval( vec_ ), E ), i );
               column( odres_ , i ) += column( expand( eval( vec_ ), E ), i );
               column( sres_  , i ) += column( expand( eval( vec_ ), E ), i );
               column( osres_ , i ) += column( expand( eval( vec_ ), E ), i );
               column( refres_, i ) += column( expand( eval( refvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) += column( expand( eval( tvec_ ), E ), i );
               column( todres_ , i ) += column( expand( eval( tvec_ ), E ), i );
               column( tsres_  , i ) += column( expand( eval( tvec_ ), E ), i );
               column( tosres_ , i ) += column( expand( eval( tvec_ ), E ), i );
               column( trefres_, i ) += column( expand( eval( trefvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Column-wise expansion with addition assignment with evaluated vector (compile time)
      {
         test_  = "Column-wise expansion with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) += column( expand<E>( eval( vec_ ) ), i );
               column( odres_ , i ) += column( expand<E>( eval( vec_ ) ), i );
               column( sres_  , i ) += column( expand<E>( eval( vec_ ) ), i );
               column( osres_ , i ) += column( expand<E>( eval( vec_ ) ), i );
               column( refres_, i ) += column( expand<E>( eval( refvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) += column( expand<E>( eval( tvec_ ) ), i );
               column( todres_ , i ) += column( expand<E>( eval( tvec_ ) ), i );
               column( tsres_  , i ) += column( expand<E>( eval( tvec_ ) ), i );
               column( tosres_ , i ) += column( expand<E>( eval( tvec_ ) ), i );
               column( trefres_, i ) += column( expand<E>( eval( trefvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Column-wise expansion with subtraction assignment
      //=====================================================================================

      // Column-wise expansion with subtraction assignment with the given vector (runtime)
      {
         test_  = "Column-wise expansion with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) -= column( expand( vec_, E ), i );
               column( odres_ , i ) -= column( expand( vec_, E ), i );
               column( sres_  , i ) -= column( expand( vec_, E ), i );
               column( osres_ , i ) -= column( expand( vec_, E ), i );
               column( refres_, i ) -= column( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) -= column( expand( tvec_, E ), i );
               column( todres_ , i ) -= column( expand( tvec_, E ), i );
               column( tsres_  , i ) -= column( expand( tvec_, E ), i );
               column( tosres_ , i ) -= column( expand( tvec_, E ), i );
               column( trefres_, i ) -= column( expand( trefvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Column-wise expansion with subtraction assignment with the given vector (compile time)
      {
         test_  = "Column-wise expansion with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) -= column( expand<E>( vec_ ), i );
               column( odres_ , i ) -= column( expand<E>( vec_ ), i );
               column( sres_  , i ) -= column( expand<E>( vec_ ), i );
               column( osres_ , i ) -= column( expand<E>( vec_ ), i );
               column( refres_, i ) -= column( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) -= column( expand<E>( tvec_ ), i );
               column( todres_ , i ) -= column( expand<E>( tvec_ ), i );
               column( tsres_  , i ) -= column( expand<E>( tvec_ ), i );
               column( tosres_ , i ) -= column( expand<E>( tvec_ ), i );
               column( trefres_, i ) -= column( expand<E>( trefvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Column-wise expansion with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Column-wise expansion with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) -= column( expand( eval( vec_ ), E ), i );
               column( odres_ , i ) -= column( expand( eval( vec_ ), E ), i );
               column( sres_  , i ) -= column( expand( eval( vec_ ), E ), i );
               column( osres_ , i ) -= column( expand( eval( vec_ ), E ), i );
               column( refres_, i ) -= column( expand( eval( refvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) -= column( expand( eval( tvec_ ), E ), i );
               column( todres_ , i ) -= column( expand( eval( tvec_ ), E ), i );
               column( tsres_  , i ) -= column( expand( eval( tvec_ ), E ), i );
               column( tosres_ , i ) -= column( expand( eval( tvec_ ), E ), i );
               column( trefres_, i ) -= column( expand( eval( trefvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Column-wise expansion with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Column-wise expansion with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) -= column( expand<E>( eval( vec_ ) ), i );
               column( odres_ , i ) -= column( expand<E>( eval( vec_ ) ), i );
               column( sres_  , i ) -= column( expand<E>( eval( vec_ ) ), i );
               column( osres_ , i ) -= column( expand<E>( eval( vec_ ) ), i );
               column( refres_, i ) -= column( expand<E>( eval( refvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) -= column( expand<E>( eval( tvec_ ) ), i );
               column( todres_ , i ) -= column( expand<E>( eval( tvec_ ) ), i );
               column( tsres_  , i ) -= column( expand<E>( eval( tvec_ ) ), i );
               column( tosres_ , i ) -= column( expand<E>( eval( tvec_ ) ), i );
               column( trefres_, i ) -= column( expand<E>( eval( trefvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Column-wise expansion with multiplication assignment
      //=====================================================================================

      // Column-wise expansion with multiplication assignment with the given vector (runtime)
      {
         test_  = "Column-wise expansion with multiplication assignment with the given vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) *= column( expand( vec_, E ), i );
               column( odres_ , i ) *= column( expand( vec_, E ), i );
               column( sres_  , i ) *= column( expand( vec_, E ), i );
               column( osres_ , i ) *= column( expand( vec_, E ), i );
               column( refres_, i ) *= column( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) *= column( expand( tvec_, E ), i );
               column( todres_ , i ) *= column( expand( tvec_, E ), i );
               column( tsres_  , i ) *= column( expand( tvec_, E ), i );
               column( tosres_ , i ) *= column( expand( tvec_, E ), i );
               column( trefres_, i ) *= column( expand( trefvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Column-wise expansion with multiplication assignment with the given vector (compile time)
      {
         test_  = "Column-wise expansion with multiplication assignment with the given vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) *= column( expand<E>( vec_ ), i );
               column( odres_ , i ) *= column( expand<E>( vec_ ), i );
               column( sres_  , i ) *= column( expand<E>( vec_ ), i );
               column( osres_ , i ) *= column( expand<E>( vec_ ), i );
               column( refres_, i ) *= column( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) *= column( expand<E>( tvec_ ), i );
               column( todres_ , i ) *= column( expand<E>( tvec_ ), i );
               column( tsres_  , i ) *= column( expand<E>( tvec_ ), i );
               column( tosres_ , i ) *= column( expand<E>( tvec_ ), i );
               column( trefres_, i ) *= column( expand<E>( trefvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Column-wise expansion with multiplication assignment with evaluated vector (runtime)
      {
         test_  = "Column-wise expansion with multiplication assignment with evaluated vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) *= column( expand( eval( vec_ ), E ), i );
               column( odres_ , i ) *= column( expand( eval( vec_ ), E ), i );
               column( sres_  , i ) *= column( expand( eval( vec_ ), E ), i );
               column( osres_ , i ) *= column( expand( eval( vec_ ), E ), i );
               column( refres_, i ) *= column( expand( eval( refvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) *= column( expand( eval( tvec_ ), E ), i );
               column( todres_ , i ) *= column( expand( eval( tvec_ ), E ), i );
               column( tsres_  , i ) *= column( expand( eval( tvec_ ), E ), i );
               column( tosres_ , i ) *= column( expand( eval( tvec_ ), E ), i );
               column( trefres_, i ) *= column( expand( eval( trefvec_ ), E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Column-wise expansion with multiplication assignment with evaluated vector (compile time)
      {
         test_  = "Column-wise expansion with multiplication assignment with evaluated vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( size_t i=0UL; i<E; ++i ) {
               column( dres_  , i ) *= column( expand<E>( eval( vec_ ) ), i );
               column( odres_ , i ) *= column( expand<E>( eval( vec_ ) ), i );
               column( sres_  , i ) *= column( expand<E>( eval( vec_ ) ), i );
               column( osres_ , i ) *= column( expand<E>( eval( vec_ ) ), i );
               column( refres_, i ) *= column( expand<E>( eval( refvec_ ) ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t i=0UL; i<vec_.size(); ++i ) {
               column( tdres_  , i ) *= column( expand<E>( eval( tvec_ ) ), i );
               column( todres_ , i ) *= column( expand<E>( eval( tvec_ ) ), i );
               column( tsres_  , i ) *= column( expand<E>( eval( tvec_ ) ), i );
               column( tosres_ , i ) *= column( expand<E>( eval( tvec_ ) ), i );
               column( trefres_, i ) *= column( expand<E>( eval( trefvec_ ) ), i );
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

      // Out-of-bounds access (invalid column index)
      {
         test_  = "Out-of-bounds column construction (invalid column index)";
         error_ = "Setup of out-of-bounds column succeeded";

         try {
            auto c = column( expand( vec_, E ), E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = column( expand<E>( vec_ ), E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = column( expand( tvec_, E ), tvec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = column( expand<E>( tvec_ ), tvec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
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
/*!\brief Skipping the column-wise dense vector expansion operation.
//
// \return void
//
// This function is called in case the column-wise dense vector expansion operation is not
// available for the given vector type \a VT.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testColumnOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the columns-wise dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the columns-wise vector expansion with plain assignment, addition assignment,
// subtraction assignment, and Schur product assignment. In case any error resulting from the
// addition or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testColumnsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION > 1 )
   {
      using blaze::expand;

      if( vec_.size() == 0UL || E == 0UL )
         return;


      std::vector<size_t> indices( E );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );

      std::vector<size_t> tindices( vec_.size() );
      std::iota( tindices.begin(), tindices.end(), 0UL );
      std::random_shuffle( tindices.begin(), tindices.end() );


      //=====================================================================================
      // Columns-wise expansion
      //=====================================================================================

      // Columns-wise expansion with the given vector (runtime)
      {
         test_  = "Columns-wise expansion with the given vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( expand( vec_, E ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( expand( vec_, E ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( expand( vec_, E ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( expand( vec_, E ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( expand( refvec_, E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) = columns( expand( tvec_, E ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) = columns( expand( tvec_, E ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) = columns( expand( tvec_, E ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) = columns( expand( tvec_, E ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) = columns( expand( trefvec_, E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Columns-wise expansion with the given vector (compile time)
      {
         test_  = "Columns-wise expansion with the given vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( expand<E>( vec_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( expand<E>( vec_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( expand<E>( vec_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( expand<E>( vec_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( expand<E>( refvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) = columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) = columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) = columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) = columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) = columns( expand<E>( trefvec_ ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Columns-wise expansion with evaluated vector (runtime)
      {
         test_  = "Columns-wise expansion with evaluated vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( expand( eval( refvec_ ), E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) = columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) = columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) = columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) = columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) = columns( expand( eval( trefvec_ ), E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Columns-wise expansion with evaluated vector (compile time)
      {
         test_  = "Columns-wise expansion with evaluated vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) = columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) = columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) = columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) = columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) = columns( expand<E>( eval( refvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) = columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) = columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) = columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) = columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) = columns( expand<E>( eval( trefvec_ ) ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Columns-wise expansion with addition assignment
      //=====================================================================================

      // Columns-wise expansion with addition assignment with the given vector (runtime)
      {
         test_  = "Columns-wise expansion with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( expand( vec_, E ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( expand( vec_, E ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( expand( vec_, E ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( expand( vec_, E ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( expand( refvec_, E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) += columns( expand( tvec_, E ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) += columns( expand( tvec_, E ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) += columns( expand( tvec_, E ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) += columns( expand( tvec_, E ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) += columns( expand( trefvec_, E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Columns-wise expansion with addition assignment with the given vector (compile time)
      {
         test_  = "Columns-wise expansion with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( expand<E>( vec_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( expand<E>( vec_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( expand<E>( vec_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( expand<E>( vec_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( expand<E>( refvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) += columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) += columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) += columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) += columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) += columns( expand<E>( trefvec_ ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Columns-wise expansion with addition assignment with evaluated vector (runtime)
      {
         test_  = "Columns-wise expansion with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( expand( eval( refvec_ ), E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) += columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) += columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) += columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) += columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) += columns( expand( eval( trefvec_ ), E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Columns-wise expansion with addition assignment with evaluated vector (compile time)
      {
         test_  = "Columns-wise expansion with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) += columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) += columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) += columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) += columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) += columns( expand<E>( eval( refvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) += columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) += columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) += columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) += columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) += columns( expand<E>( eval( trefvec_ ) ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Columns-wise expansion with subtraction assignment
      //=====================================================================================

      // Columns-wise expansion with subtraction assignment with the given vector (runtime)
      {
         test_  = "Columns-wise expansion with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( expand( vec_, E ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( expand( vec_, E ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( expand( vec_, E ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( expand( vec_, E ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( expand( refvec_, E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) -= columns( expand( tvec_, E ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) -= columns( expand( tvec_, E ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) -= columns( expand( tvec_, E ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) -= columns( expand( tvec_, E ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) -= columns( expand( trefvec_, E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Columns-wise expansion with subtraction assignment with the given vector (compile time)
      {
         test_  = "Columns-wise expansion with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( expand<E>( vec_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( expand<E>( vec_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( expand<E>( vec_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( expand<E>( vec_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( expand<E>( refvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) -= columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) -= columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) -= columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) -= columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) -= columns( expand<E>( trefvec_ ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Columns-wise expansion with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Columns-wise expansion with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( expand( eval( refvec_ ), E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) -= columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) -= columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) -= columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) -= columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) -= columns( expand( eval( trefvec_ ), E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Columns-wise expansion with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Columns-wise expansion with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) -= columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) -= columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) -= columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) -= columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) -= columns( expand<E>( eval( refvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) -= columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) -= columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) -= columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) -= columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) -= columns( expand<E>( eval( trefvec_ ) ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Columns-wise expansion with Schur product assignment
      //=====================================================================================

      // Columns-wise expansion with Schur product assignment with the given vector (runtime)
      {
         test_  = "Columns-wise expansion with Schur product assignment with the given vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( expand( vec_, E ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( expand( vec_, E ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( expand( vec_, E ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( expand( vec_, E ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( expand( refvec_, E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) %= columns( expand( tvec_, E ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) %= columns( expand( tvec_, E ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) %= columns( expand( tvec_, E ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) %= columns( expand( tvec_, E ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) %= columns( expand( trefvec_, E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Columns-wise expansion with Schur product assignment with the given vector (compile time)
      {
         test_  = "Columns-wise expansion with Schur product assignment with the given vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( expand<E>( vec_ ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( expand<E>( vec_ ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( expand<E>( vec_ ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( expand<E>( vec_ ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( expand<E>( refvec_ ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) %= columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) %= columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) %= columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) %= columns( expand<E>( tvec_ ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) %= columns( expand<E>( trefvec_ ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Columns-wise expansion with Schur product assignment with evaluated vector (runtime)
      {
         test_  = "Columns-wise expansion with Schur product assignment with evaluated vector (runtime)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( expand( eval( vec_ ), E ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( expand( eval( refvec_ ), E ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) %= columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) %= columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) %= columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) %= columns( expand( eval( tvec_ ), E ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) %= columns( expand( eval( trefvec_ ), E ), &tindices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Columns-wise expansion with Schur product assignment with evaluated vector (compile time)
      {
         test_  = "Columns-wise expansion with Schur product assignment with evaluated vector (compile time)";
         error_ = "Failed Schur product assignment";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               columns( dres_  , &indices[index], n ) %= columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( odres_ , &indices[index], n ) %= columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( sres_  , &indices[index], n ) %= columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( osres_ , &indices[index], n ) %= columns( expand<E>( eval( vec_ ) ), &indices[index], n );
               columns( refres_, &indices[index], n ) %= columns( expand<E>( eval( refvec_ ) ), &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( size_t index=0UL, n=0UL; index<tindices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, tindices.size() - index );
               columns( tdres_  , &tindices[index], n ) %= columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( todres_ , &tindices[index], n ) %= columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( tsres_  , &tindices[index], n ) %= columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( tosres_ , &tindices[index], n ) %= columns( expand<E>( eval( tvec_ ) ), &tindices[index], n );
               columns( trefres_, &tindices[index], n ) %= columns( expand<E>( eval( trefvec_ ) ), &tindices[index], n );
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

      // Out-of-bounds access (invalid column index; initializer_list)
      {
         test_  = "Out-of-bounds column selection construction (invalid column index; initializer_list)";
         error_ = "Setup of out-of-bounds column selection succeeded";

         try {
            auto c = columns( expand( vec_, E ), { E } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( expand<E>( vec_ ), { E } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( expand( tvec_, E ), { tvec_.size() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( expand<E>( tvec_ ), { tvec_.size() } );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
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
            auto c = columns( expand( vec_, E ), [index=E]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( expand<E>( vec_ ), [index=E]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( expand( tvec_, E ), [index=tvec_.size()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
                << "   Result:\n" << c << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid column access index" );
         }

         try {
            auto c = columns( expand<E>( tvec_ ), [index=tvec_.size()]( size_t ){ return index; }, 1UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
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
/*!\brief Skipping the columns-wise dense vector expansion operation.
//
// \return void
//
// This function is called in case the columns-wise dense vector expansion operation is not
// available for the given vector type \a VT.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testColumnsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the band-wise dense vector expansion operation.
//
// \return void
// \exception std::runtime_error Expansion error detected.
//
// This function tests the band-wise vector expansion with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment. In case any error resulting from the
// addition or the subsequent assignment is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testBandOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_BAND_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BAND_OPERATION > 1 )
   {
      using blaze::expand;

      if( vec_.size() == 0UL || E == 0UL )
         return;


      //=====================================================================================
      // Band-wise expansion
      //=====================================================================================

      // Band-wise expansion with the given vector (runtime)
      {
         test_  = "Band-wise expansion with the given vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) = band( expand( vec_, E )   , i );
               band( odres_ , i ) = band( expand( vec_, E )   , i );
               band( sres_  , i ) = band( expand( vec_, E )   , i );
               band( osres_ , i ) = band( expand( vec_, E )   , i );
               band( refres_, i ) = band( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) = band( expand( tvec_, E )   , j );
               band( todres_ , j ) = band( expand( tvec_, E )   , j );
               band( tsres_  , j ) = band( expand( tvec_, E )   , j );
               band( tosres_ , j ) = band( expand( tvec_, E )   , j );
               band( trefres_, j ) = band( expand( trefvec_, E ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Band-wise expansion with the given vector (compile time)
      {
         test_  = "Band-wise expansion with the given vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) = band( expand<E>( vec_ )   , i );
               band( odres_ , i ) = band( expand<E>( vec_ )   , i );
               band( sres_  , i ) = band( expand<E>( vec_ )   , i );
               band( osres_ , i ) = band( expand<E>( vec_ )   , i );
               band( refres_, i ) = band( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) = band( expand<E>( tvec_ )   , j );
               band( todres_ , j ) = band( expand<E>( tvec_ )   , j );
               band( tsres_  , j ) = band( expand<E>( tvec_ )   , j );
               band( tosres_ , j ) = band( expand<E>( tvec_ )   , j );
               band( trefres_, j ) = band( expand<E>( trefvec_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Band-wise expansion with evaluated vector (runtime)
      {
         test_  = "Band-wise expansion with evaluated vector (runtime)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) = band( expand( vec_, E )   , i );
               band( odres_ , i ) = band( expand( vec_, E )   , i );
               band( sres_  , i ) = band( expand( vec_, E )   , i );
               band( osres_ , i ) = band( expand( vec_, E )   , i );
               band( refres_, i ) = band( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) = band( expand( tvec_, E )   , j );
               band( todres_ , j ) = band( expand( tvec_, E )   , j );
               band( tsres_  , j ) = band( expand( tvec_, E )   , j );
               band( tosres_ , j ) = band( expand( tvec_, E )   , j );
               band( trefres_, j ) = band( expand( trefvec_, E ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Band-wise expansion with evaluated vector (compile time)
      {
         test_  = "Band-wise expansion with evaluated vector (compile time)";
         error_ = "Failed expansion operation";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) = band( expand<E>( vec_ )   , i );
               band( odres_ , i ) = band( expand<E>( vec_ )   , i );
               band( sres_  , i ) = band( expand<E>( vec_ )   , i );
               band( osres_ , i ) = band( expand<E>( vec_ )   , i );
               band( refres_, i ) = band( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) = band( expand<E>( tvec_ )   , j );
               band( todres_ , j ) = band( expand<E>( tvec_ )   , j );
               band( tsres_  , j ) = band( expand<E>( tvec_ )   , j );
               band( tosres_ , j ) = band( expand<E>( tvec_ )   , j );
               band( trefres_, j ) = band( expand<E>( trefvec_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Band-wise expansion with addition assignment
      //=====================================================================================

      // Band-wise expansion with addition assignment with the given vector (runtime)
      {
         test_  = "Band-wise expansion with addition assignment with the given vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) += band( expand( vec_, E )   , i );
               band( odres_ , i ) += band( expand( vec_, E )   , i );
               band( sres_  , i ) += band( expand( vec_, E )   , i );
               band( osres_ , i ) += band( expand( vec_, E )   , i );
               band( refres_, i ) += band( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) += band( expand( tvec_, E )   , j );
               band( todres_ , j ) += band( expand( tvec_, E )   , j );
               band( tsres_  , j ) += band( expand( tvec_, E )   , j );
               band( tosres_ , j ) += band( expand( tvec_, E )   , j );
               band( trefres_, j ) += band( expand( trefvec_, E ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Band-wise expansion with addition assignment with the given vector (compile time)
      {
         test_  = "Band-wise expansion with addition assignment with the given vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) += band( expand<E>( vec_ )   , i );
               band( odres_ , i ) += band( expand<E>( vec_ )   , i );
               band( sres_  , i ) += band( expand<E>( vec_ )   , i );
               band( osres_ , i ) += band( expand<E>( vec_ )   , i );
               band( refres_, i ) += band( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) += band( expand<E>( tvec_ )   , j );
               band( todres_ , j ) += band( expand<E>( tvec_ )   , j );
               band( tsres_  , j ) += band( expand<E>( tvec_ )   , j );
               band( tosres_ , j ) += band( expand<E>( tvec_ )   , j );
               band( trefres_, j ) += band( expand<E>( trefvec_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Band-wise expansion with addition assignment with evaluated vector (runtime)
      {
         test_  = "Band-wise expansion with addition assignment with evaluated vector (runtime)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) += band( expand( vec_, E )   , i );
               band( odres_ , i ) += band( expand( vec_, E )   , i );
               band( sres_  , i ) += band( expand( vec_, E )   , i );
               band( osres_ , i ) += band( expand( vec_, E )   , i );
               band( refres_, i ) += band( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) += band( expand( tvec_, E )   , j );
               band( todres_ , j ) += band( expand( tvec_, E )   , j );
               band( tsres_  , j ) += band( expand( tvec_, E )   , j );
               band( tosres_ , j ) += band( expand( tvec_, E )   , j );
               band( trefres_, j ) += band( expand( trefvec_, E ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Band-wise expansion with addition assignment with evaluated vector (compile time)
      {
         test_  = "Band-wise expansion with addition assignment with evaluated vector (compile time)";
         error_ = "Failed addition assignment";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) += band( expand<E>( vec_ )   , i );
               band( odres_ , i ) += band( expand<E>( vec_ )   , i );
               band( sres_  , i ) += band( expand<E>( vec_ )   , i );
               band( osres_ , i ) += band( expand<E>( vec_ )   , i );
               band( refres_, i ) += band( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) += band( expand<E>( tvec_ )   , j );
               band( todres_ , j ) += band( expand<E>( tvec_ )   , j );
               band( tsres_  , j ) += band( expand<E>( tvec_ )   , j );
               band( tosres_ , j ) += band( expand<E>( tvec_ )   , j );
               band( trefres_, j ) += band( expand<E>( trefvec_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Band-wise expansion with subtraction assignment
      //=====================================================================================

      // Band-wise expansion with subtraction assignment with the given vector (runtime)
      {
         test_  = "Band-wise expansion with subtraction assignment with the given vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) -= band( expand( vec_, E )   , i );
               band( odres_ , i ) -= band( expand( vec_, E )   , i );
               band( sres_  , i ) -= band( expand( vec_, E )   , i );
               band( osres_ , i ) -= band( expand( vec_, E )   , i );
               band( refres_, i ) -= band( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) -= band( expand( tvec_, E )   , j );
               band( todres_ , j ) -= band( expand( tvec_, E )   , j );
               band( tsres_  , j ) -= band( expand( tvec_, E )   , j );
               band( tosres_ , j ) -= band( expand( tvec_, E )   , j );
               band( trefres_, j ) -= band( expand( trefvec_, E ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Band-wise expansion with subtraction assignment with the given vector (compile time)
      {
         test_  = "Band-wise expansion with subtraction assignment with the given vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) -= band( expand<E>( vec_ )   , i );
               band( odres_ , i ) -= band( expand<E>( vec_ )   , i );
               band( sres_  , i ) -= band( expand<E>( vec_ )   , i );
               band( osres_ , i ) -= band( expand<E>( vec_ )   , i );
               band( refres_, i ) -= band( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) -= band( expand<E>( tvec_ )   , j );
               band( todres_ , j ) -= band( expand<E>( tvec_ )   , j );
               band( tsres_  , j ) -= band( expand<E>( tvec_ )   , j );
               band( tosres_ , j ) -= band( expand<E>( tvec_ )   , j );
               band( trefres_, j ) -= band( expand<E>( trefvec_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Band-wise expansion with subtraction assignment with evaluated vector (runtime)
      {
         test_  = "Band-wise expansion with subtraction assignment with evaluated vector (runtime)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) -= band( expand( vec_, E )   , i );
               band( odres_ , i ) -= band( expand( vec_, E )   , i );
               band( sres_  , i ) -= band( expand( vec_, E )   , i );
               band( osres_ , i ) -= band( expand( vec_, E )   , i );
               band( refres_, i ) -= band( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) -= band( expand( tvec_, E )   , j );
               band( todres_ , j ) -= band( expand( tvec_, E )   , j );
               band( tsres_  , j ) -= band( expand( tvec_, E )   , j );
               band( tosres_ , j ) -= band( expand( tvec_, E )   , j );
               band( trefres_, j ) -= band( expand( trefvec_, E ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Band-wise expansion with subtraction assignment with evaluated vector (compile time)
      {
         test_  = "Band-wise expansion with subtraction assignment with evaluated vector (compile time)";
         error_ = "Failed subtraction assignment";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) -= band( expand<E>( vec_ )   , i );
               band( odres_ , i ) -= band( expand<E>( vec_ )   , i );
               band( sres_  , i ) -= band( expand<E>( vec_ )   , i );
               band( osres_ , i ) -= band( expand<E>( vec_ )   , i );
               band( refres_, i ) -= band( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) -= band( expand<E>( tvec_ )   , j );
               band( todres_ , j ) -= band( expand<E>( tvec_ )   , j );
               band( tsres_  , j ) -= band( expand<E>( tvec_ )   , j );
               band( tosres_ , j ) -= band( expand<E>( tvec_ )   , j );
               band( trefres_, j ) -= band( expand<E>( trefvec_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }


      //=====================================================================================
      // Band-wise expansion with multiplication assignment
      //=====================================================================================

      // Band-wise expansion with multiplication assignment with the given vector (runtime)
      {
         test_  = "Band-wise expansion with multiplication assignment with the given vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) *= band( expand( vec_, E )   , i );
               band( odres_ , i ) *= band( expand( vec_, E )   , i );
               band( sres_  , i ) *= band( expand( vec_, E )   , i );
               band( osres_ , i ) *= band( expand( vec_, E )   , i );
               band( refres_, i ) *= band( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) *= band( expand( tvec_, E )   , j );
               band( todres_ , j ) *= band( expand( tvec_, E )   , j );
               band( tsres_  , j ) *= band( expand( tvec_, E )   , j );
               band( tosres_ , j ) *= band( expand( tvec_, E )   , j );
               band( trefres_, j ) *= band( expand( trefvec_, E ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Band-wise expansion with multiplication assignment with the given vector (compile time)
      {
         test_  = "Band-wise expansion with multiplication assignment with the given vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) *= band( expand<E>( vec_ )   , i );
               band( odres_ , i ) *= band( expand<E>( vec_ )   , i );
               band( sres_  , i ) *= band( expand<E>( vec_ )   , i );
               band( osres_ , i ) *= band( expand<E>( vec_ )   , i );
               band( refres_, i ) *= band( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) *= band( expand<E>( tvec_ )   , j );
               band( todres_ , j ) *= band( expand<E>( tvec_ )   , j );
               band( tsres_  , j ) *= band( expand<E>( tvec_ )   , j );
               band( tosres_ , j ) *= band( expand<E>( tvec_ )   , j );
               band( trefres_, j ) *= band( expand<E>( trefvec_ ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Band-wise expansion with multiplication assignment with evaluated vector (runtime)
      {
         test_  = "Band-wise expansion with multiplication assignment with evaluated vector (runtime)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) *= band( expand( vec_, E )   , i );
               band( odres_ , i ) *= band( expand( vec_, E )   , i );
               band( sres_  , i ) *= band( expand( vec_, E )   , i );
               band( osres_ , i ) *= band( expand( vec_, E )   , i );
               band( refres_, i ) *= band( expand( refvec_, E ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) *= band( expand( tvec_, E )   , j );
               band( todres_ , j ) *= band( expand( tvec_, E )   , j );
               band( tsres_  , j ) *= band( expand( tvec_, E )   , j );
               band( tosres_ , j ) *= band( expand( tvec_, E )   , j );
               band( trefres_, j ) *= band( expand( trefvec_, E ), j );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT>( ex );
         }

         checkTransposeResults<TVT>();
      }

      // Band-wise expansion with multiplication assignment with evaluated vector (compile time)
      {
         test_  = "Band-wise expansion with multiplication assignment with evaluated vector (compile time)";
         error_ = "Failed multiplication assignment";

         try {
            initResults();
            for( ptrdiff_t i=1UL-vec_.size(); i<E; ++i ) {
               band( dres_  , i ) *= band( expand<E>( vec_ )   , i );
               band( odres_ , i ) *= band( expand<E>( vec_ )   , i );
               band( sres_  , i ) *= band( expand<E>( vec_ )   , i );
               band( osres_ , i ) *= band( expand<E>( vec_ )   , i );
               band( refres_, i ) *= band( expand<E>( refvec_ ), i );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT>( ex );
         }

         checkResults<VT>();

         try {
            initTransposeResults();
            for( ptrdiff_t j=1UL-E; j<tvec_.size(); ++j ) {
               band( tdres_  , j ) *= band( expand<E>( tvec_ )   , j );
               band( todres_ , j ) *= band( expand<E>( tvec_ )   , j );
               band( tsres_  , j ) *= band( expand<E>( tvec_ )   , j );
               band( tosres_ , j ) *= band( expand<E>( tvec_ )   , j );
               band( trefres_, j ) *= band( expand<E>( trefvec_ ), j );
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

      // Out-of-bounds access (invalid lower band index)
      {
         test_  = "Out-of-bounds band construction (invalid lower band index)";
         error_ = "Setup of out-of-bounds band succeeded";

         try {
            auto b = band( expand( vec_, E ), -vec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( expand<E>( vec_ ), -vec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( expand( tvec_, E ), -E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( expand<E>( tvec_ ), -E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
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
            auto b = band( expand( vec_, E ), E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( expand<E>( vec_ ), E );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( VT ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( expand( tvec_, E ), tvec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
                << "   Result:\n" << b << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ex ) {
            checkExceptionMessage( ex, "Invalid band access index" );
         }

         try {
            auto b = band( expand<E>( tvec_ ), tvec_.size() );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: " << error_ << "\n"
                << " Details:\n"
                << "   Random seed = " << blaze::getSeed() << "\n"
                << "   Dense vector type:\n"
                << "     " << typeid( TVT ).name() << "\n"
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
/*!\brief Skipping the band-wise dense vector expansion operation.
//
// \return void
//
// This function is called in case the band-wise dense vector expansion operation is not available
// for the given vector type \a VT.
*/
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::testBandOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the customized dense vector expansion operation.
//
// \param op The custom operation to be tested.
// \param name The human-readable name of the operation.
// \return void
// \exception std::runtime_error Operation error detected.
//
// This function tests the vector expansion operation with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment in combination with
// a custom operation. In case any error resulting from the expansion or the subsequent assignment
// is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT    // Type of the dense vector
        , size_t E >     // Compile time expansion
template< typename OP >  // Type of the custom operation
void OperationTest<VT,E>::testCustomOperation( OP op, const std::string& name )
{
   using blaze::expand;


   //=====================================================================================
   // Customized expansion operation
   //=====================================================================================

   // Customized expansion operation with the given vector (runtime)
   {
      test_  = "Customized expansion operation with the given vector (runtime)";
      error_ = "Failed expansion operation";

      try {
         initResults();
         dres_   = op( expand( vec_, E ) );
         odres_  = op( expand( vec_, E ) );
         sres_   = op( expand( vec_, E ) );
         osres_  = op( expand( vec_, E ) );
         refres_ = op( expand( refvec_, E ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   = op( expand( tvec_, E ) );
         todres_  = op( expand( tvec_, E ) );
         tsres_   = op( expand( tvec_, E ) );
         tosres_  = op( expand( tvec_, E ) );
         trefres_ = op( expand( trefvec_, E ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized expansion operation with the given vector (compile time)
   {
      test_  = "Customized expansion operation with the given vector (compile time)";
      error_ = "Failed expansion operation";

      try {
         initResults();
         dres_   = op( expand<E>( vec_ ) );
         odres_  = op( expand<E>( vec_ ) );
         sres_   = op( expand<E>( vec_ ) );
         osres_  = op( expand<E>( vec_ ) );
         refres_ = op( expand<E>( refvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   = op( expand<E>( tvec_ ) );
         todres_  = op( expand<E>( tvec_ ) );
         tsres_   = op( expand<E>( tvec_ ) );
         tosres_  = op( expand<E>( tvec_ ) );
         trefres_ = op( expand<E>( trefvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized expansion operation with evaluated vector (runtime)
   {
      test_  = "Customized expansion operation with evaluated vector (runtime)";
      error_ = "Failed expansion operation";

      try {
         initResults();
         dres_   = op( expand( eval( vec_ ), E ) );
         odres_  = op( expand( eval( vec_ ), E ) );
         sres_   = op( expand( eval( vec_ ), E ) );
         osres_  = op( expand( eval( vec_ ), E ) );
         refres_ = op( expand( eval( refvec_ ), E ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   = op( expand( eval( tvec_ ), E ) );
         todres_  = op( expand( eval( tvec_ ), E ) );
         tsres_   = op( expand( eval( tvec_ ), E ) );
         tosres_  = op( expand( eval( tvec_ ), E ) );
         trefres_ = op( expand( eval( trefvec_ ), E ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized expansion operation with evaluated vector (compile time)
   {
      test_  = "Customized expansion operation with evaluated vector (compile time)";
      error_ = "Failed expansion operation";

      try {
         initResults();
         dres_   = op( expand<E>( eval( vec_ ) ) );
         odres_  = op( expand<E>( eval( vec_ ) ) );
         sres_   = op( expand<E>( eval( vec_ ) ) );
         osres_  = op( expand<E>( eval( vec_ ) ) );
         refres_ = op( expand<E>( eval( refvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   = op( expand<E>( eval( tvec_ ) ) );
         todres_  = op( expand<E>( eval( tvec_ ) ) );
         tsres_   = op( expand<E>( eval( tvec_ ) ) );
         tosres_  = op( expand<E>( eval( tvec_ ) ) );
         trefres_ = op( expand<E>( eval( trefvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }


   //=====================================================================================
   // Customized expansion with addition assignment
   //=====================================================================================

   // Customized expansion with addition assignment with the given vector (runtime)
   {
      test_  = "Customized expansion with addition assignment with the given vector (runtime)";
      error_ = "Failed addition assignment";

      try {
         initResults();
         dres_   += op( expand( vec_, E ) );
         odres_  += op( expand( vec_, E ) );
         sres_   += op( expand( vec_, E ) );
         osres_  += op( expand( vec_, E ) );
         refres_ += op( expand( refvec_, E ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   += op( expand( tvec_, E ) );
         todres_  += op( expand( tvec_, E ) );
         tsres_   += op( expand( tvec_, E ) );
         tosres_  += op( expand( tvec_, E ) );
         trefres_ += op( expand( trefvec_, E ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized expansion with addition assignment with the given vector (compile time)
   {
      test_  = "Customized expansion with addition assignment with the given vector (compile time)";
      error_ = "Failed addition assignment";

      try {
         initResults();
         dres_   += op( expand<E>( vec_ ) );
         odres_  += op( expand<E>( vec_ ) );
         sres_   += op( expand<E>( vec_ ) );
         osres_  += op( expand<E>( vec_ ) );
         refres_ += op( expand<E>( refvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   += op( expand<E>( tvec_ ) );
         todres_  += op( expand<E>( tvec_ ) );
         tsres_   += op( expand<E>( tvec_ ) );
         tosres_  += op( expand<E>( tvec_ ) );
         trefres_ += op( expand<E>( trefvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized expansion with addition assignment with evaluated vector (runtime)
   {
      test_  = "Customized expansion with addition assignment with evaluated vector (runtime)";
      error_ = "Failed addition assignment";

      try {
         initResults();
         dres_   += op( expand( eval( vec_ ), E ) );
         odres_  += op( expand( eval( vec_ ), E ) );
         sres_   += op( expand( eval( vec_ ), E ) );
         osres_  += op( expand( eval( vec_ ), E ) );
         refres_ += op( expand( eval( refvec_ ), E ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   += op( expand( eval( tvec_ ), E ) );
         todres_  += op( expand( eval( tvec_ ), E ) );
         tsres_   += op( expand( eval( tvec_ ), E ) );
         tosres_  += op( expand( eval( tvec_ ), E ) );
         trefres_ += op( expand( eval( trefvec_ ), E ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized expansion with addition assignment with evaluated vector (compile time)
   {
      test_  = "Customized expansion with addition assignment with evaluated vector (compile time)";
      error_ = "Failed addition assignment";

      try {
         initResults();
         dres_   += op( expand<E>( eval( vec_ ) ) );
         odres_  += op( expand<E>( eval( vec_ ) ) );
         sres_   += op( expand<E>( eval( vec_ ) ) );
         osres_  += op( expand<E>( eval( vec_ ) ) );
         refres_ += op( expand<E>( eval( refvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   += op( expand<E>( eval( tvec_ ) ) );
         todres_  += op( expand<E>( eval( tvec_ ) ) );
         tsres_   += op( expand<E>( eval( tvec_ ) ) );
         tosres_  += op( expand<E>( eval( tvec_ ) ) );
         trefres_ += op( expand<E>( eval( trefvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }


   //=====================================================================================
   // Customized expansion with subtraction assignment
   //=====================================================================================

   // Customized expansion with subtraction assignment with the given vector (runtime)
   {
      test_  = "Customized expansion with subtraction assignment with the given vector (runtime)";
      error_ = "Failed subtraction assignment";

      try {
         initResults();
         dres_   -= op( expand( vec_, E ) );
         odres_  -= op( expand( vec_, E ) );
         sres_   -= op( expand( vec_, E ) );
         osres_  -= op( expand( vec_, E ) );
         refres_ -= op( expand( refvec_, E ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   -= op( expand( tvec_, E ) );
         todres_  -= op( expand( tvec_, E ) );
         tsres_   -= op( expand( tvec_, E ) );
         tosres_  -= op( expand( tvec_, E ) );
         trefres_ -= op( expand( trefvec_, E ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized expansion with subtraction assignment with the given vector (compile time)
   {
      test_  = "Customized expansion with subtraction assignment with the given vector (compile time)";
      error_ = "Failed subtraction assignment";

      try {
         initResults();
         dres_   -= op( expand<E>( vec_ ) );
         odres_  -= op( expand<E>( vec_ ) );
         sres_   -= op( expand<E>( vec_ ) );
         osres_  -= op( expand<E>( vec_ ) );
         refres_ -= op( expand<E>( refvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   -= op( expand<E>( tvec_ ) );
         todres_  -= op( expand<E>( tvec_ ) );
         tsres_   -= op( expand<E>( tvec_ ) );
         tosres_  -= op( expand<E>( tvec_ ) );
         trefres_ -= op( expand<E>( trefvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized expansion with subtraction assignment with evaluated vector (runtime)
   {
      test_  = "Customized expansion with subtraction assignment with evaluated vector (runtime)";
      error_ = "Failed subtraction assignment";

      try {
         initResults();
         dres_   -= op( expand( eval( vec_ ), E ) );
         odres_  -= op( expand( eval( vec_ ), E ) );
         sres_   -= op( expand( eval( vec_ ), E ) );
         osres_  -= op( expand( eval( vec_ ), E ) );
         refres_ -= op( expand( eval( refvec_ ), E ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   -= op( expand( eval( tvec_ ), E ) );
         todres_  -= op( expand( eval( tvec_ ), E ) );
         tsres_   -= op( expand( eval( tvec_ ), E ) );
         tosres_  -= op( expand( eval( tvec_ ), E ) );
         trefres_ -= op( expand( eval( trefvec_ ), E ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized expansion with subtraction assignment with evaluated vector (compile time)
   {
      test_  = "Customized expansion with subtraction assignment with evaluated vector (compile time)";
      error_ = "Failed subtraction assignment";

      try {
         initResults();
         dres_   -= op( expand<E>( eval( vec_ ) ) );
         odres_  -= op( expand<E>( eval( vec_ ) ) );
         sres_   -= op( expand<E>( eval( vec_ ) ) );
         osres_  -= op( expand<E>( eval( vec_ ) ) );
         refres_ -= op( expand<E>( eval( refvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   -= op( expand<E>( eval( tvec_ ) ) );
         todres_  -= op( expand<E>( eval( tvec_ ) ) );
         tsres_   -= op( expand<E>( eval( tvec_ ) ) );
         tosres_  -= op( expand<E>( eval( tvec_ ) ) );
         trefres_ -= op( expand<E>( eval( trefvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }


   //=====================================================================================
   // Customized expansion with Schur product assignment
   //=====================================================================================

   // Customized expansion with Schur product assignment with the given vector (runtime)
   {
      test_  = "Customized expansion with Schur product assignment with the given vector (runtime)";
      error_ = "Failed Schur product assignment";

      try {
         initResults();
         dres_   %= op( expand( vec_, E ) );
         odres_  %= op( expand( vec_, E ) );
         sres_   %= op( expand( vec_, E ) );
         osres_  %= op( expand( vec_, E ) );
         refres_ %= op( expand( refvec_, E ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   %= op( expand( tvec_, E ) );
         todres_  %= op( expand( tvec_, E ) );
         tsres_   %= op( expand( tvec_, E ) );
         tosres_  %= op( expand( tvec_, E ) );
         trefres_ %= op( expand( trefvec_, E ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized expansion with Schur product assignment with the given vector (compile time)
   {
      test_  = "Customized expansion with Schur product assignment with the given vector (compile time)";
      error_ = "Failed Schur product assignment";

      try {
         initResults();
         dres_   %= op( expand<E>( vec_ ) );
         odres_  %= op( expand<E>( vec_ ) );
         sres_   %= op( expand<E>( vec_ ) );
         osres_  %= op( expand<E>( vec_ ) );
         refres_ %= op( expand<E>( refvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   %= op( expand<E>( tvec_ ) );
         todres_  %= op( expand<E>( tvec_ ) );
         tsres_   %= op( expand<E>( tvec_ ) );
         tosres_  %= op( expand<E>( tvec_ ) );
         trefres_ %= op( expand<E>( trefvec_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized expansion with Schur product assignment with evaluated vector (runtime)
   {
      test_  = "Customized expansion with Schur product assignment with evaluated vector (runtime)";
      error_ = "Failed Schur product assignment";

      try {
         initResults();
         dres_   %= op( expand( eval( vec_ ), E ) );
         odres_  %= op( expand( eval( vec_ ), E ) );
         sres_   %= op( expand( eval( vec_ ), E ) );
         osres_  %= op( expand( eval( vec_ ), E ) );
         refres_ %= op( expand( eval( refvec_ ), E ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   %= op( expand( eval( tvec_ ), E ) );
         todres_  %= op( expand( eval( tvec_ ), E ) );
         tsres_   %= op( expand( eval( tvec_ ), E ) );
         tosres_  %= op( expand( eval( tvec_ ), E ) );
         trefres_ %= op( expand( eval( trefvec_ ), E ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT>( ex );
      }

      checkTransposeResults<TVT>();
   }

   // Customized expansion with Schur product assignment with evaluated vector (compile time)
   {
      test_  = "Customized expansion with Schur product assignment with evaluated vector (compile time)";
      error_ = "Failed Schur product assignment";

      try {
         initResults();
         dres_   %= op( expand<E>( eval( vec_ ) ) );
         odres_  %= op( expand<E>( eval( vec_ ) ) );
         sres_   %= op( expand<E>( eval( vec_ ) ) );
         osres_  %= op( expand<E>( eval( vec_ ) ) );
         refres_ %= op( expand<E>( eval( refvec_ ) ) );
      }
      catch( std::exception& ex ) {
         convertException<VT>( ex );
      }

      checkResults<VT>();

      try {
         initTransposeResults();
         tdres_   %= op( expand<E>( eval( tvec_ ) ) );
         todres_  %= op( expand<E>( eval( tvec_ ) ) );
         tsres_   %= op( expand<E>( eval( tvec_ ) ) );
         tosres_  %= op( expand<E>( eval( tvec_ ) ) );
         trefres_ %= op( expand<E>( eval( trefvec_ ) ) );
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
// template argument \a T indicates the type of the vector operand used for the computations.
*/
template< typename VT   // Type of the dense vector
        , size_t E >    // Compile time expansion
template< typename T >  // Type of the vector operand
void OperationTest<VT,E>::checkResults()
{
   using blaze::IsRowVector;

   if( !isEqual( dres_, refres_ ) || !isEqual( odres_, refres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result matrix detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense " << ( IsRowVector<T>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
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
          << "   Dense " << ( IsRowVector<T>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
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
// results. The template argument \a T indicates the type of the vector operand used for the
// computations.
*/
template< typename VT   // Type of the dense vector
        , size_t E >    // Compile time expansion
template< typename T >  // Type of the vector operand
void OperationTest<VT,E>::checkTransposeResults()
{
   using blaze::IsRowVector;

   if( !isEqual( tdres_, trefres_ ) || !isEqual( todres_, trefres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect dense result matrix detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense " << ( IsRowVector<T>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( T ).name() << "\n"
          << "   Transpose result:\n" << tdres_ << "\n"
          << "   Transpose result with opposite storage order:\n" << todres_ << "\n"
          << "   Expected result:\n" << trefres_ << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( !isEqual( tsres_, trefres_ ) || !isEqual( tosres_, trefres_ ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect sparse result matrix detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Dense " << ( IsRowVector<T>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
          << "     " << typeid( T ).name() << "\n"
          << "   Transpose result:\n" << tsres_ << "\n"
          << "   Transpose result with opposite storage order:\n" << tosres_ << "\n"
          << "   Expected result:\n" << trefres_ << "\n";
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
template< typename VT   // Type of the dense vector
        , size_t E >    // Compile time expansion
void OperationTest<VT,E>::checkExceptionMessage( const std::exception& ex, const std::string& message )
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
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::initResults()
{
   const blaze::UnderlyingBuiltin_t<DRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<DRE> max( randmax );

   resize( dres_, size( vec_ ), E );
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
template< typename VT  // Type of the dense vector
        , size_t E >   // Compile time expansion
void OperationTest<VT,E>::initTransposeResults()
{
   const blaze::UnderlyingBuiltin_t<DRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<DRE> max( randmax );

   resize( tdres_, E, size( tvec_ ) );
   randomize( tdres_, min, max );

   todres_  = tdres_;
   tsres_   = tdres_;
   tosres_  = tdres_;
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
// test. The template argument \a T indicates the type of vector operand used for the computations.
*/
template< typename VT   // Type of the dense vector
        , size_t E >    // Compile time expansion
template< typename T >  // Type of the vector operand
void OperationTest<VT,E>::convertException( const std::exception& ex )
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
/*!\brief Testing the expansion operation for a specific vector type.
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
/*!\brief Macro for the definition of a dense vector expansion operation test case.
*/
#define DEFINE_DVECEXPAND_OPERATION_TEST( VT ) \
   extern template class blazetest::mathtest::operations::dvecexpand::OperationTest<VT>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a dense vector expansion operation test case.
*/
#define RUN_DVECEXPAND_OPERATION_TEST( C ) \
   blazetest::mathtest::operations::dvecexpand::runTest( C )
/*! \endcond */
//*************************************************************************************************

} // namespace dvecexpand

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
