//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/dvecdvecmin/OperationTest.h
//  \brief Header file for the dense vector/dense vector minimum operation test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_DVECDVECMIN_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_DVECDVECMIN_OPERATIONTEST_H_


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
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/Functors.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDivisor.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/traits/MapTrait.h>
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

namespace dvecdvecmin {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the dense vector/dense vector minimum operation test.
//
// This class template represents one particular vector minimum test between two vectors of
// a particular type. The two template arguments \a VT1 and \a VT2 represent the types of the
// left-hand side and right-hand side vector, respectively.
*/
template< typename VT1    // Type of the left-hand side dense vector
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

   using DRE  = blaze::MapTrait_t<VT1,VT2,blaze::Min>;    //!< Dense result type
   using TDRE = blaze::MapTrait_t<TVT1,TVT2,blaze::Min>;  //!< Transpose dense result type
   using DET  = blaze::ElementType_t<DRE>;                //!< Element type of the dense result

   using SRE  = blaze::CompressedVector<DET,TF>;  //!< Sparse result type
   using TSRE = blaze::TransposeType_t<SRE>;      //!< Transpose sparse result type
   using SET  = blaze::ElementType_t<SRE>;        //!< Element type of the sparse result

   using RT  = blaze::DynamicVector<blaze::ElementType_t<DRE>,TF>;  //!< Reference type
   using TRT = blaze::TransposeType_t<RT>;                          //!< Transpose reference type

   //! Type of the vector/vector minimum expression
   using VecVecMinExprType =
      blaze::RemoveCVRef_t< decltype( min( std::declval<VT1>(), std::declval<VT2>() ) ) >;

   //! Type of the transpose vector/transpose vector minimum expression
   using TVecTVecMinExprType =
      blaze::RemoveCVRef_t< decltype( min( std::declval<TVT1>(), std::declval<TVT2>() ) ) >;
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
   VT1  lhs_;      //!< The left-hand side dense vector.
   VT2  rhs_;      //!< The right-hand side dense vector.
   DRE  dres_;     //!< The dense vector for the result of the vector minimum.
   SRE  sres_;     //!< The sparse vector for the result of the vector minimum.
   RT   ref_;      //!< The reference vector.
   RT   refres_;   //!< The reference result.
   TVT1 tlhs_;     //!< The transpose left-hand side vector.
   TVT2 trhs_;     //!< The transpose right-hand side vector.
   TDRE tdres_;    //!< The dense vector for the result of the transpose vector minimum.
   TSRE tsres_;    //!< The sparse vector for the result of the transpose vector minimum.
   TRT  tref_;     //!< The transpose reference vector.
   TRT  trefres_;  //!< The transpose reference result.

   std::string test_;   //!< Label of the currently performed test.
   std::string error_;  //!< Description of the current error type.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT1  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT2  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TVT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TVT2 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( RT   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TRT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( DRE  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( SRE  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( TDRE );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( TSRE );

   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , VT2  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , RT   );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TVT1, TVT2 );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( TVT1, TRT  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , DRE  );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG( VT1 , SRE  );
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
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RT , blaze::TransposeType_t<TRT>  );

   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG     ( VecVecMinExprType, blaze::ResultType_t<VecVecMinExprType>    );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_DIFFERENT_TRANSPOSE_FLAG( VecVecMinExprType, blaze::TransposeType_t<VecVecMinExprType> );

   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_SAME_TRANSPOSE_FLAG     ( TVecTVecMinExprType, blaze::ResultType_t<TVecTVecMinExprType>    );
   BLAZE_CONSTRAINT_VECTORS_MUST_HAVE_DIFFERENT_TRANSPOSE_FLAG( TVecTVecMinExprType, blaze::TransposeType_t<TVecTVecMinExprType> );

   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( VecVecMinExprType  , blaze::BaseType_t<VecVecMinExprType  > );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( TVecTVecMinExprType, blaze::BaseType_t<TVecTVecMinExprType> );
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
/*!\brief Constructor for the dense vector/dense vector minimum operation test.
//
// \param creator1 The creator for the left-hand side dense vector of the vector minimum.
// \param creator2 The creator for the right-hand side dense vector of the vector minimum.
// \exception std::runtime_error Operation error detected.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
OperationTest<VT1,VT2>::OperationTest( const Creator<VT1>& creator1, const Creator<VT2>& creator2 )
   : lhs_( creator1() )     // The left-hand side dense vector
   , rhs_( creator2() )     // The right-hand side dense vector
   , dres_()                // The dense vector for the result of the vector minimum
   , sres_()                // The sparse vector for the result of the vector minimum
   , ref_()                 // The reference vector
   , refres_()              // The reference result
   , tlhs_( trans(lhs_) )   // The transpose left-hand side vector
   , trhs_( trans(rhs_) )   // The transpose right-hand side vector
   , tdres_()               // The dense vector for the result of the transpose vector minimum
   , tsres_()               // The sparse vector for the result of the transpose vector minimum
   , tref_()                // The transpose reference vector
   , trefres_()             // The transpose reference result
   , test_()                // Label of the currently performed test
   , error_()               // Description of the current error type
{
   using namespace blaze;

   using Scalar = UnderlyingScalar_t<DET>;

   if( lhs_.size() != rhs_.size() ) {
      throw std::runtime_error( "Non-matching operands detected" );
   }

   ref_.resize( lhs_.size() );
   tref_.resize( tlhs_.size() );
   for( size_t i=0UL; i<lhs_.size(); ++i ) {
      ref_ [i] = min( lhs_ [i], rhs_ [i] );
      tref_[i] = min( tlhs_[i], trhs_[i] );
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


   //=====================================================================================
   // Testing the evaluation with the given vectors
   //=====================================================================================

   try
   {
      const auto res   ( evaluate( min( lhs_, rhs_ ) ) );
      const auto refres( evaluate( ref_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the given vectors\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense " << ( IsRowVector<VT1>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
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
   catch( std::exception& ex ) {
      convertException<VT1,VT2>( ex );
   }

   try
   {
      const auto res   ( evaluate( min( eval( lhs_ ), eval( rhs_ ) ) ) );
      const auto refres( evaluate( eval( ref_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated vectors\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense " << ( IsRowVector<VT1>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
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
   catch( std::exception& ex ) {
      convertException<VT1,VT2>( ex );
   }


   //=====================================================================================
   // Testing the evaluation with the transpose types
   //=====================================================================================

   try
   {
      const auto res   ( evaluate( min( tlhs_, trhs_ ) ) );
      const auto refres( evaluate( tref_ ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with the transpose vectors\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense " << ( IsRowVector<TVT1>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
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
   catch( std::exception& ex ) {
      convertException<TVT1,TVT2>( ex );
   }

   try
   {
      const auto res   ( evaluate( min( eval( tlhs_ ), eval( trhs_ ) ) ) );
      const auto refres( evaluate( eval( tref_ ) ) );

      if( !isEqual( res, refres ) ) {
         std::ostringstream oss;
         oss << " Test: Evaluation with evaluated transpose vectors\n"
             << " Error: Failed evaluation\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense " << ( IsRowVector<TVT1>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
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
   catch( std::exception& ex ) {
      convertException<TVT1,TVT2>( ex );
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
template< typename VT1    // Type of the left-hand side dense vector
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

      if( !equal( min( lhs_, rhs_ )[n], ref_[n] ) ||
          !equal( min( lhs_, rhs_ ).at(n), ref_.at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of minimum expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( min( lhs_, eval( rhs_ ) )[n], ref_[n] ) ||
          !equal( min( lhs_, eval( rhs_ ) ).at(n), ref_.at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated minimum expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( min( eval( lhs_ ), rhs_ )[n], ref_[n] ) ||
          !equal( min( eval( lhs_ ), rhs_ ).at(n), ref_.at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated minimum expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( min( eval( lhs_ ), eval( rhs_ ) )[n], ref_[n] ) |
          !equal( min( eval( lhs_ ), eval( rhs_ ) ).at(n), ref_.at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated minimum expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
             << "     " << typeid( VT1 ).name() << "\n"
             << "   Right-hand side dense vector type:\n"
             << "     " << typeid( VT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      min( lhs_, rhs_ ).at( lhs_.size() );

      std::ostringstream oss;
      oss << " Test : Checked element access of minimum expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Left-hand side dense vector type:\n"
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

      if( !equal( min( tlhs_, trhs_ )[n], tref_[n] ) ||
          !equal( min( tlhs_, trhs_ ).at(n), tref_.at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of transpose minimum expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Transpose left-hand side dense vector type:\n"
             << "     " << typeid( TVT1 ).name() << "\n"
             << "   Transpose right-hand side dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( min( tlhs_, eval( trhs_ ) )[n], tref_[n] ) ||
          !equal( min( tlhs_, eval( trhs_ ) ).at(n), tref_.at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of right evaluated transpose minimum expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Transpose left-hand side dense vector type:\n"
             << "     " << typeid( TVT1 ).name() << "\n"
             << "   Transpose right-hand side dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( min( eval( tlhs_ ), trhs_ )[n], tref_[n] ) ||
          !equal( min( eval( tlhs_ ), trhs_ ).at(n), tref_.at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of left evaluated transpose minimum expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Transpose left-hand side dense vector type:\n"
             << "     " << typeid( TVT1 ).name() << "\n"
             << "   Transpose right-hand side dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( !equal( min( eval( tlhs_ ), eval( trhs_ ) )[n], tref_[n] ) ||
          !equal( min( eval( tlhs_ ), eval( trhs_ ) ).at(n), tref_.at(n) ) ) {
         std::ostringstream oss;
         oss << " Test : Element access of fully evaluated transpose minimum expression\n"
             << " Error: Unequal resulting elements at index " << n << " detected\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Transpose left-hand side dense vector type:\n"
             << "     " << typeid( TVT1 ).name() << "\n"
             << "   Transpose right-hand side dense vector type:\n"
             << "     " << typeid( TVT2 ).name() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      min( tlhs_, trhs_ ).at( tlhs_.size() );

      std::ostringstream oss;
      oss << " Test : Checked element access of transpose minimum expression\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   Transpose left-hand side dense vector type:\n"
          << "     " << typeid( TVT1 ).name() << "\n"
          << "   Transpose right-hand side dense vector type:\n"
          << "     " << typeid( TVT2 ).name() << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the plain dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the vector minimum with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case
// any error resulting from the minimum operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testBasicOperation()
{
#if BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
   if( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 1 )
   {
      //=====================================================================================
      // Minimum
      //=====================================================================================

      // Minimum with the given vectors
      {
         test_  = "Minimum with the given vectors";
         error_ = "Failed minimum operation";

         try {
            initResults();
            dres_   = min( lhs_, rhs_ );
            sres_   = min( lhs_, rhs_ );
            refres_ = ref_;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = min( tlhs_, trhs_ );
            tsres_   = min( tlhs_, trhs_ );
            trefres_ = tref_;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Minimum with evaluated vectors
      {
         test_  = "Minimum with evaluated vectors";
         error_ = "Failed minimum operation";

         try {
            initResults();
            dres_   = min( eval( lhs_ ), eval( rhs_ ) );
            sres_   = min( eval( lhs_ ), eval( rhs_ ) );
            refres_ = eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = min( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   = min( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ = eval( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Minimum with addition assignment
      //=====================================================================================

      // Minimum with addition assignment with the given vectors
      {
         test_  = "Minimum with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += min( lhs_, rhs_ );
            sres_   += min( lhs_, rhs_ );
            refres_ += ref_;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += min( tlhs_, trhs_ );
            tsres_   += min( tlhs_, trhs_ );
            trefres_ += tref_;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Minimum with addition assignment with the given vectors
      {
         test_  = "Minimum with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += min( eval( lhs_ ), eval( rhs_ ) );
            sres_   += min( eval( lhs_ ), eval( rhs_ ) );
            refres_ += eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += min( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   += min( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ += eval( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Minimum with subtraction assignment
      //=====================================================================================

      // Minimum with subtraction assignment with the given vectors
      {
         test_  = "Minimum with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= min( lhs_, rhs_ );
            sres_   -= min( lhs_, rhs_ );
            refres_ -= ref_;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= min( tlhs_, trhs_ );
            tsres_   -= min( tlhs_, trhs_ );
            trefres_ -= tref_;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Minimum with subtraction assignment with evaluated vectors
      {
         test_  = "Minimum with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= min( eval( lhs_ ), eval( rhs_ ) );
            sres_   -= min( eval( lhs_ ), eval( rhs_ ) );
            refres_ -= eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= min( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   -= min( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ -= eval( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Minimum with multiplication assignment
      //=====================================================================================

      // Minimum with multiplication assignment with the given vectors
      {
         test_  = "Minimum with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= min( lhs_, rhs_ );
            sres_   *= min( lhs_, rhs_ );
            refres_ *= ref_;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= min( tlhs_, trhs_ );
            tsres_   *= min( tlhs_, trhs_ );
            trefres_ *= tref_;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Minimum with multiplication assignment with evaluated vectors
      {
         test_  = "Minimum with multiplication assignment with evaluated vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= min( eval( lhs_ ), eval( rhs_ ) );
            sres_   *= min( eval( lhs_ ), eval( rhs_ ) );
            refres_ *= eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= min( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   *= min( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ *= eval( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Minimum with division assignment
      //=====================================================================================

      if( blaze::isDivisor( min( lhs_, rhs_ ) ) )
      {
         // Minimum with division assignment with the given vectors
         {
            test_  = "Minimum with division assignment with the given vectors";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= min( lhs_, rhs_ );
               sres_   /= min( lhs_, rhs_ );
               refres_ /= ref_;
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkResults<VT1,VT2>();

            try {
               initTransposeResults();
               tdres_   /= min( tlhs_, trhs_ );
               tsres_   /= min( tlhs_, trhs_ );
               trefres_ /= tref_;
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkTransposeResults<TVT1,TVT2>();
         }

         // Minimum with division assignment with evaluated vectors
         {
            test_  = "Minimum with division assignment with evaluated vectors";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= min( eval( lhs_ ), eval( rhs_ ) );
               sres_   /= min( eval( lhs_ ), eval( rhs_ ) );
               refres_ /= eval( ref_ );
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkResults<VT1,VT2>();

            try {
               initTransposeResults();
               tdres_   /= min( eval( tlhs_ ), eval( trhs_ ) );
               tsres_   /= min( eval( tlhs_ ), eval( trhs_ ) );
               trefres_ /= eval( tref_ );
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkTransposeResults<TVT1,TVT2>();
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the negated dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the negated vector minimum with plain assignment, addition assignment,
// subtraction assignment, and multiplication assignment, and division assignment. In case
// any error resulting from the minimum operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testNegatedOperation()
{
#if BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
   if( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 1 )
   {
      //=====================================================================================
      // Negated minimum
      //=====================================================================================

      // Negated minimum with the given vectors
      {
         test_  = "Negated minimum with the given types";
         error_ = "Failed minimum operation";

         try {
            initResults();
            dres_   = -min( lhs_, rhs_ );
            sres_   = -min( lhs_, rhs_ );
            refres_ = -ref_;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = -min( tlhs_, trhs_ );
            tsres_   = -min( tlhs_, trhs_ );
            trefres_ = -tref_;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Negated minimum with evaluated vectors
      {
         test_  = "Negated minimum with evaluated vectors";
         error_ = "Failed minimum operation";

         try {
            initResults();
            dres_   = -min( eval( lhs_ ), eval( rhs_ ) );
            sres_   = -min( eval( lhs_ ), eval( rhs_ ) );
            refres_ = -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = -min( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   = -min( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ = -eval( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Negated minimum with addition assignment
      //=====================================================================================

      // Negated minimum with addition assignment with the given vectors
      {
         test_  = "Negated minimum with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -min( lhs_, rhs_ );
            sres_   += -min( lhs_, rhs_ );
            refres_ += -ref_;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += -min( tlhs_, trhs_ );
            tsres_   += -min( tlhs_, trhs_ );
            trefres_ += -tref_;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Negated minimum with addition assignment with evaluated vectors
      {
         test_  = "Negated minimum with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += -min( eval( lhs_ ), eval( rhs_ ) );
            sres_   += -min( eval( lhs_ ), eval( rhs_ ) );
            refres_ += -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += -min( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   += -min( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ += -eval( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Negated minimum with subtraction assignment
      //=====================================================================================

      // Negated minimum with subtraction assignment with the given vectors
      {
         test_  = "Negated minimum with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -min( lhs_, rhs_ );
            sres_   -= -min( lhs_, rhs_ );
            refres_ -= -ref_;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();
         try {
            initTransposeResults();
            tdres_   -= -min( tlhs_, trhs_ );
            tsres_   -= -min( tlhs_, trhs_ );
            trefres_ -= -tref_;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Negated minimum with subtraction assignment with evaluated vectors
      {
         test_  = "Negated minimum with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= -min( eval( lhs_ ), eval( rhs_ ) );
            sres_   -= -min( eval( lhs_ ), eval( rhs_ ) );
            refres_ -= -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= -min( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   -= -min( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ -= -eval( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Negated minimum with multiplication assignment
      //=====================================================================================

      // Negated minimum with multiplication assignment with the given vectors
      {
         test_  = "Negated minimum with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= -min( lhs_, rhs_ );
            sres_   *= -min( lhs_, rhs_ );
            refres_ *= -ref_;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= -min( tlhs_, trhs_ );
            tsres_   *= -min( tlhs_, trhs_ );
            trefres_ *= -tref_;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Negated minimum with multiplication assignment with evaluated vectors
      {
         test_  = "Negated minimum with multiplication assignment with evaluated vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= -min( eval( lhs_ ), eval( rhs_ ) );
            sres_   *= -min( eval( lhs_ ), eval( rhs_ ) );
            refres_ *= -eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= -min( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   *= -min( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ *= -eval( tref_);
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Negated minimum with division assignment
      //=====================================================================================

      if( blaze::isDivisor( min( lhs_, rhs_ ) ) )
      {
         // Negated minimum with division assignment with the given vectors
         {
            test_  = "Negated minimum with division assignment with the given vectors";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= -min( lhs_, rhs_ );
               sres_   /= -min( lhs_, rhs_ );
               refres_ /= -ref_;
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkResults<VT1,VT2>();

            try {
               initTransposeResults();
               tdres_   /= -min( tlhs_, trhs_ );
               tsres_   /= -min( tlhs_, trhs_ );
               trefres_ /= -tref_;
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkTransposeResults<TVT1,TVT2>();
         }

         // Negated minimum with division assignment with evaluated vectors
         {
            test_  = "Negated minimum with division assignment with evaluated vectors";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= -min( eval( lhs_ ), eval( rhs_ ) );
               sres_   /= -min( eval( lhs_ ), eval( rhs_ ) );
               refres_ /= -eval( ref_ );
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkResults<VT1,VT2>();

            try {
               initTransposeResults();
               tdres_   /= -min( eval( tlhs_ ), eval( trhs_ ) );
               tsres_   /= -min( eval( tlhs_ ), eval( trhs_ ) );
               trefres_ /= -eval( tref_ );
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkTransposeResults<TVT1,TVT2>();
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the scaled dense vector/dense vector minimum operation.
//
// \param scalar The scalar value.
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the scaled vector minimum with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case any
// error resulting from the minimum operation or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
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
      // Self-scaling (v*=s)
      //=====================================================================================

      // Self-scaling (v*=s)
      {
         test_ = "Self-scaling (v*=s)";

         try {
            dres_   = min( lhs_, rhs_ );
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
            dres_   = min( lhs_, rhs_ );
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
            dres_   = min( lhs_, rhs_ );
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
            dres_   = min( lhs_, rhs_ );
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
            dres_   = min( lhs_, rhs_ );
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
      // Scaled minimum (s*OP)
      //=====================================================================================

      // Scaled minimum with the given vectors
      {
         test_  = "Scaled minimum with the given vectors (s*OP)";
         error_ = "Failed minimum operation";

         try {
            initResults();
            dres_   = scalar * min( lhs_, rhs_ );
            sres_   = scalar * min( lhs_, rhs_ );
            refres_ = scalar * ref_;
         }
         catch( std::exception& ex ) {
             convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = scalar * min( tlhs_, trhs_ );
            tsres_   = scalar * min( tlhs_, trhs_ );
            trefres_ = scalar * tref_;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled minimum with evaluated vectors
      {
         test_ = "Scaled minimum with evaluated vectors (s*OP)";
         error_ = "Failed minimum operation";

         try {
            initResults();
            dres_   = scalar * min( eval( lhs_ ), eval( rhs_ ) );
            sres_   = scalar * min( eval( lhs_ ), eval( rhs_ ) );
            refres_ = scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = scalar * min( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   = scalar * min( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ = scalar * eval( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled minimum (OP*s)
      //=====================================================================================

      // Scaled minimum with the given vectors
      {
         test_  = "Scaled minimum with the given vectors (OP*s)";
         error_ = "Failed minimum operation";

         try {
            initResults();
            dres_   = min( lhs_, rhs_ ) * scalar;
            sres_   = min( lhs_, rhs_ ) * scalar;
            refres_ = ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = min( tlhs_, trhs_ ) * scalar;
            tsres_   = min( tlhs_, trhs_ ) * scalar;
            trefres_ = tref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled minimum with evaluated vectors
      {
         test_  = "Scaled minimum with evaluated vectors (OP*s)";
         error_ = "Failed minimum operation";

         try {
            initResults();
            dres_   = min( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            sres_   = min( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            refres_ = eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = min( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            tsres_   = min( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            trefres_ = eval( tref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled minimum (OP/s)
      //=====================================================================================

      // Scaled minimum with the given vectors
      {
         test_  = "Scaled minimum with the given vectors (OP/s)";
         error_ = "Failed minimum operation";

         try {
            initResults();
            dres_   = min( lhs_, rhs_ ) / scalar;
            sres_   = min( lhs_, rhs_ ) / scalar;
            refres_ = ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = min( tlhs_, trhs_ ) / scalar;
            tsres_   = min( tlhs_, trhs_ ) / scalar;
            trefres_ = tref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled minimum with evaluated vectors
      {
         test_  = "Scaled minimum with evaluated vectors (OP/s)";
         error_ = "Failed minimum operation";

         try {
            initResults();
            dres_   = min( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            sres_   = min( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            refres_ = eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   = min( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            tsres_   = min( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            trefres_ = eval( tref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled minimum with addition assignment (s*OP)
      //=====================================================================================

      // Scaled minimum with addition assignment with the given vectors
      {
         test_  = "Scaled minimum with addition assignment with the given vectors (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * min( lhs_, rhs_ );
            sres_   += scalar * min( lhs_, rhs_ );
            refres_ += scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += scalar * min( tlhs_, trhs_ );
            tsres_   += scalar * min( tlhs_, trhs_ );
            trefres_ += scalar * tref_;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled minimum with addition assignment with evaluated vectors
      {
         test_  = "Scaled minimum with addition assignment with evaluated vectors (s*OP)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += scalar * min( eval( lhs_ ), eval( rhs_ ) );
            sres_   += scalar * min( eval( lhs_ ), eval( rhs_ ) );
            refres_ += scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += scalar * min( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   += scalar * min( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ += scalar * eval( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled minimum with addition assignment (OP*s)
      //=====================================================================================

      // Scaled minimum with addition assignment with the given vectors
      {
         test_  = "Scaled minimum with addition assignment with the given vectors (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += min( lhs_, rhs_ ) * scalar;
            sres_   += min( lhs_, rhs_ ) * scalar;
            refres_ += ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += min( tlhs_, trhs_ ) * scalar;
            tsres_   += min( tlhs_, trhs_ ) * scalar;
            trefres_ += tref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled minimum with addition assignment with evaluated vectors
      {
         test_  = "Scaled minimum with addition assignment with evaluated vectors (OP*s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += min( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            sres_   += min( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            refres_ += eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += min( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            tsres_   += min( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            trefres_ += eval( tref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled minimum with addition assignment (OP/s)
      //=====================================================================================

      // Scaled minimum with addition assignment with the given vectors
      {
         test_  = "Scaled minimum with addition assignment with the given vectors (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += min( lhs_, rhs_ ) / scalar;
            sres_   += min( lhs_, rhs_ ) / scalar;
            refres_ += ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += min( tlhs_, trhs_ ) / scalar;
            tsres_   += min( tlhs_, trhs_ ) / scalar;
            trefres_ += tref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled minimum with addition assignment with evaluated vectors
      {
         test_  = "Scaled minimum with addition assignment with evaluated vectors (OP/s)";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            dres_   += min( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            sres_   += min( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            refres_ += eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   += min( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            tsres_   += min( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            trefres_ += eval( tref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled minimum with subtraction assignment (s*OP)
      //=====================================================================================

      // Scaled minimum with subtraction assignment with the given vectors
      {
         test_  = "Scaled minimum with subtraction assignment with the given vectors (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * min( lhs_, rhs_ );
            sres_   -= scalar * min( lhs_, rhs_ );
            refres_ -= scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= scalar * min( tlhs_, trhs_ );
            tsres_   -= scalar * min( tlhs_, trhs_ );
            trefres_ -= scalar * tref_;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled minimum with subtraction assignment with evaluated vectors
      {
         test_  = "Scaled minimum with subtraction assignment with evaluated vectors (s*OP)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= scalar * min( eval( lhs_ ), eval( rhs_ ) );
            sres_   -= scalar * min( eval( lhs_ ), eval( rhs_ ) );
            refres_ -= scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= scalar * min( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   -= scalar * min( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ -= scalar * eval( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled minimum with subtraction assignment (OP*s)
      //=====================================================================================

      // Scaled minimum with subtraction assignment with the given vectors
      {
         test_  = "Scaled minimum with subtraction assignment with the given vectors (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= min( lhs_, rhs_ ) * scalar;
            sres_   -= min( lhs_, rhs_ ) * scalar;
            refres_ -= ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= min( tlhs_, trhs_ ) * scalar;
            tsres_   -= min( tlhs_, trhs_ ) * scalar;
            trefres_ -= tref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled minimum with subtraction assignment with evaluated vectors
      {
         test_  = "Scaled minimum with subtraction assignment with evaluated vectors (OP*s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= min( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            sres_   -= min( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            refres_ -= eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= min( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            tsres_   -= min( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            trefres_ -= eval( tref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled minimum with subtraction assignment (OP/s)
      //=====================================================================================

      // Scaled minimum with subtraction assignment with the given vectors
      {
         test_  = "Scaled minimum with subtraction assignment with the given vectors (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= min( lhs_, rhs_ ) / scalar;
            sres_   -= min( lhs_, rhs_ ) / scalar;
            refres_ -= ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= min( tlhs_, trhs_ ) / scalar;
            tsres_   -= min( tlhs_, trhs_ ) / scalar;
            trefres_ -= tref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled minimum with subtraction assignment with evaluated vectors
      {
         test_  = "Scaled minimum with subtraction assignment with evaluated vectors (OP/s)";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            dres_   -= min( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            sres_   -= min( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            refres_ -= eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   -= min( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            tsres_   -= min( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            trefres_ -= eval( tref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled minimum with multiplication assignment (s*OP)
      //=====================================================================================

      // Scaled minimum with multiplication assignment with the given vectors
      {
         test_  = "Scaled minimum with multiplication assignment with the given vectors (s*OP)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= scalar * min( lhs_, rhs_ );
            sres_   *= scalar * min( lhs_, rhs_ );
            refres_ *= scalar * ref_;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= scalar * min( tlhs_, trhs_ );
            tsres_   *= scalar * min( tlhs_, trhs_ );
            trefres_ *= scalar * tref_;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled minimum with multiplication assignment with evaluated vectors
      {
         test_  = "Scaled minimum with multiplication assignment with evaluated vectors (s*OP)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= scalar * min( eval( lhs_ ), eval( rhs_ ) );
            sres_   *= scalar * min( eval( lhs_ ), eval( rhs_ ) );
            refres_ *= scalar * eval( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= scalar * min( eval( tlhs_ ), eval( trhs_ ) );
            tsres_   *= scalar * min( eval( tlhs_ ), eval( trhs_ ) );
            trefres_ *= scalar * eval( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled minimum with multiplication assignment (OP*s)
      //=====================================================================================

      // Scaled minimum with multiplication assignment with the given vectors
      {
         test_  = "Scaled minimum with multiplication assignment with the given vectors (OP*s)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= min( lhs_, rhs_ ) * scalar;
            sres_   *= min( lhs_, rhs_ ) * scalar;
            refres_ *= ref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= min( tlhs_, trhs_ ) * scalar;
            tsres_   *= min( tlhs_, trhs_ ) * scalar;
            trefres_ *= tref_ * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled minimum with multiplication assignment with evaluated vectors
      {
         test_  = "Scaled minimum with multiplication assignment with evaluated vectors (OP*s)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= min( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            sres_   *= min( eval( lhs_ ), eval( rhs_ ) ) * scalar;
            refres_ *= eval( ref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= min( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            tsres_   *= min( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
            trefres_ *= eval( tref_ ) * scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled minimum with multiplication assignment (OP/s)
      //=====================================================================================

      // Scaled minimum with multiplication assignment with the given vectors
      {
         test_  = "Scaled minimum with multiplication assignment with the given vectors (OP/s)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= min( lhs_, rhs_ ) / scalar;
            sres_   *= min( lhs_, rhs_ ) / scalar;
            refres_ *= ref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= min( tlhs_, trhs_ ) / scalar;
            tsres_   *= min( tlhs_, trhs_ ) / scalar;
            trefres_ *= tref_ / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Scaled minimum with multiplication assignment with evaluated vectors
      {
         test_  = "Scaled minimum with multiplication assignment with evaluated vectors (OP/s)";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            dres_   *= min( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            sres_   *= min( eval( lhs_ ), eval( rhs_ ) ) / scalar;
            refres_ *= eval( ref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   *= min( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            tsres_   *= min( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
            trefres_ *= eval( tref_ ) / scalar;
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Scaled minimum with division assignment (s*OP)
      //=====================================================================================

      if( blaze::isDivisor( min( lhs_, rhs_ ) ) )
      {
         // Scaled minimum with division assignment with the given vectors
         {
            test_  = "Scaled minimum with division assignment with the given vectors (s*OP)";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= scalar * min( lhs_, rhs_ );
               sres_   /= scalar * min( lhs_, rhs_ );
               refres_ /= scalar * ref_;
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkResults<VT1,VT2>();

            try {
               initTransposeResults();
               tdres_   /= scalar * min( tlhs_, trhs_ );
               tsres_   /= scalar * min( tlhs_, trhs_ );
               trefres_ /= scalar * tref_;
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkTransposeResults<TVT1,TVT2>();
         }

         // Scaled minimum with division assignment with evaluated vectors
         {
            test_  = "Scaled minimum with division assignment with evaluated vectors (s*OP)";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= scalar * min( eval( lhs_ ), eval( rhs_ ) );
               sres_   /= scalar * min( eval( lhs_ ), eval( rhs_ ) );
               refres_ /= scalar * eval( ref_ );
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkResults<VT1,VT2>();

            try {
               initTransposeResults();
               tdres_   /= scalar * min( eval( tlhs_ ), eval( trhs_ ) );
               tsres_   /= scalar * min( eval( tlhs_ ), eval( trhs_ ) );
               trefres_ /= scalar * eval( tref_ );
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkTransposeResults<TVT1,TVT2>();
         }
      }


      //=====================================================================================
      // Scaled minimum with division assignment (OP*s)
      //=====================================================================================

      if( blaze::isDivisor( min( lhs_, rhs_ ) ) )
      {
         // Scaled minimum with division assignment with the given vectors
         {
            test_  = "Scaled minimum with division assignment with the given vectors (OP*s)";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= min( lhs_, rhs_ ) * scalar;
               sres_   /= min( lhs_, rhs_ ) * scalar;
               refres_ /= ref_ * scalar;
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkResults<VT1,VT2>();

            try {
               initTransposeResults();
               tdres_   /= min( tlhs_, trhs_ ) * scalar;
               tsres_   /= min( tlhs_, trhs_ ) * scalar;
               trefres_ /= tref_ * scalar;
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkTransposeResults<TVT1,TVT2>();
         }

         // Scaled minimum with division assignment with evaluated vectors
         {
            test_  = "Scaled minimum with division assignment with evaluated vectors (OP*s)";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= min( eval( lhs_ ), eval( rhs_ ) ) * scalar;
               sres_   /= min( eval( lhs_ ), eval( rhs_ ) ) * scalar;
               refres_ /= eval( ref_ ) * scalar;
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkResults<VT1,VT2>();

            try {
               initTransposeResults();
               tdres_   /= min( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
               tsres_   /= min( eval( tlhs_ ), eval( trhs_ ) ) * scalar;
               trefres_ /= eval( tref_ ) * scalar;
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkTransposeResults<TVT1,TVT2>();
         }
      }


      //=====================================================================================
      // Scaled minimum with division assignment (OP/s)
      //=====================================================================================

      if( blaze::isDivisor( ( min( lhs_, rhs_ ) ) / scalar ) )
      {
         // Scaled minimum with division assignment with the given vectors
         {
            test_  = "Scaled minimum with division assignment with the given vectors (OP/s)";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= min( lhs_, rhs_ ) / scalar;
               sres_   /= min( lhs_, rhs_ ) / scalar;
               refres_ /= ref_ / scalar;
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkResults<VT1,VT2>();

            try {
               initTransposeResults();
               tdres_   /= min( tlhs_, trhs_ ) / scalar;
               tsres_   /= min( tlhs_, trhs_ ) / scalar;
               trefres_ /= tref_ / scalar;
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkTransposeResults<TVT1,TVT2>();
         }

         // Scaled minimum with division assignment with evaluated vectors
         {
            test_  = "Scaled minimum with division assignment with evaluated vectors (OP/s)";
            error_ = "Failed division assignment operation";

            try {
               initResults();
               dres_   /= min( eval( lhs_ ), eval( rhs_ ) ) / scalar;
               sres_   /= min( eval( lhs_ ), eval( rhs_ ) ) / scalar;
               refres_ /= eval( ref_ ) / scalar;
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkResults<VT1,VT2>();

            try {
               initTransposeResults();
               tdres_   /= min( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
               tsres_   /= min( eval( tlhs_ ), eval( trhs_ ) ) / scalar;
               trefres_ /= eval( tref_ ) / scalar;
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkTransposeResults<TVT1,TVT2>();
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the transpose dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the transpose vector minimum with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case any
// error resulting from the minimum operation or the subsequent assignment is detected, a
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
      // Transpose minimum
      //=====================================================================================

      // Transpose minimum with the given vectors
      {
         test_  = "Transpose minimum with the given vectors";
         error_ = "Failed minimum operation";

         try {
            initTransposeResults();
            tdres_   = trans( min( lhs_, rhs_ ) );
            tsres_   = trans( min( lhs_, rhs_ ) );
            trefres_ = trans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   = trans( min( tlhs_, trhs_ ) );
            sres_   = trans( min( tlhs_, trhs_ ) );
            refres_ = trans( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Transpose minimum with evaluated vectors
      {
         test_  = "Transpose minimum with evaluated vectors";
         error_ = "Failed minimum operation";

         try {
            initTransposeResults();
            tdres_   = trans( min( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   = trans( min( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ = trans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   = trans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   = trans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ = trans( eval( tref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Transpose minimum with addition assignment
      //=====================================================================================

      // Transpose minimum with addition assignment with the given vectors
      {
         test_  = "Transpose minimum with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initTransposeResults();
            tdres_   += trans( min( lhs_, rhs_ ) );
            tsres_   += trans( min( lhs_, rhs_ ) );
            trefres_ += trans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   += trans( min( tlhs_, trhs_ ) );
            sres_   += trans( min( tlhs_, trhs_ ) );
            refres_ += trans( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Transpose minimum with addition assignment with evaluated vectors
      {
         test_  = "Transpose minimum with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initTransposeResults();
            tdres_   += trans( min( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   += trans( min( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ += trans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   += trans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   += trans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ += trans( eval( tref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Transpose minimum with subtraction assignment
      //=====================================================================================

      // Transpose minimum with subtraction assignment with the given vectors
      {
         test_  = "Transpose minimum with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initTransposeResults();
            tdres_   -= trans( min( lhs_, rhs_ ) );
            tsres_   -= trans( min( lhs_, rhs_ ) );
            trefres_ -= trans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   -= trans( min( tlhs_, trhs_ ) );
            sres_   -= trans( min( tlhs_, trhs_ ) );
            refres_ -= trans( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Transpose minimum with subtraction assignment with evaluated vectors
      {
         test_  = "Transpose minimum with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initTransposeResults();
            tdres_   -= trans( min( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   -= trans( min( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ -= trans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   -= trans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   -= trans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ -= trans( eval( tref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Transpose minimum with multiplication assignment
      //=====================================================================================

      // Transpose minimum with multiplication assignment with the given vectors
      {
         test_  = "Transpose minimum with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initTransposeResults();
            tdres_   *= trans( min( lhs_, rhs_ ) );
            tsres_   *= trans( min( lhs_, rhs_ ) );
            trefres_ *= trans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   *= trans( min( tlhs_, trhs_ ) );
            sres_   *= trans( min( tlhs_, trhs_ ) );
            refres_ *= trans( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Transpose minimum with multiplication assignment with evaluated vectors
      {
         test_  = "Transpose minimum with multiplication assignment with evaluated vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initTransposeResults();
            tdres_   *= trans( min( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   *= trans( min( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ *= trans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   *= trans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   *= trans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ *= trans( eval( tref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Transpose minimum with division assignment
      //=====================================================================================

      if( blaze::isDivisor( min( lhs_, rhs_ ) ) )
      {
         // Transpose minimum with division assignment with the given vectors
         {
            test_  = "Transpose minimum with division assignment with the given vectors";
            error_ = "Failed division assignment operation";

            try {
               initTransposeResults();
               tdres_   /= trans( min( lhs_, rhs_ ) );
               tsres_   /= trans( min( lhs_, rhs_ ) );
               trefres_ /= trans( ref_ );
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkTransposeResults<VT1,VT2>();

            try {
               initResults();
               dres_   /= trans( min( tlhs_, trhs_ ) );
               sres_   /= trans( min( tlhs_, trhs_ ) );
               refres_ /= trans( tref_ );
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkResults<TVT1,TVT2>();
         }

         // Transpose minimum with division assignment with evaluated vectors
         {
            test_  = "Transpose minimum with division assignment with evaluated vectors";
            error_ = "Failed division assignment operation";

            try {
               initTransposeResults();
               tdres_   /= trans( min( eval( lhs_ ), eval( rhs_ ) ) );
               tsres_   /= trans( min( eval( lhs_ ), eval( rhs_ ) ) );
               trefres_ /= trans( eval( ref_ ) );
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkTransposeResults<VT1,VT2>();

            try {
               initResults();
               dres_   /= trans( min( eval( tlhs_ ), eval( trhs_ ) ) );
               sres_   /= trans( min( eval( tlhs_ ), eval( trhs_ ) ) );
               refres_ /= trans( eval( tref_ ) );
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkResults<TVT1,TVT2>();
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the conjugate transpose dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the conjugate transpose vector minimum with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment. In
// case any error resulting from the minimum operation or the subsequent assignment is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testCTransOperation()
{
#if BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION > 1 )
   {
      //=====================================================================================
      // Conjugate transpose minimum
      //=====================================================================================

      // Conjugate transpose minimum with the given vectors
      {
         test_  = "Conjugate transpose minimum with the given vectors";
         error_ = "Failed minimum operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( min( lhs_, rhs_ ) );
            tsres_   = ctrans( min( lhs_, rhs_ ) );
            trefres_ = ctrans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   = ctrans( min( tlhs_, trhs_ ) );
            sres_   = ctrans( min( tlhs_, trhs_ ) );
            refres_ = ctrans( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Conjugate transpose minimum with evaluated vectors
      {
         test_  = "Conjugate transpose minimum with evaluated vectors";
         error_ = "Failed minimum operation";

         try {
            initTransposeResults();
            tdres_   = ctrans( min( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   = ctrans( min( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ = ctrans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   = ctrans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   = ctrans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ = ctrans( eval( tref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Conjugate transpose minimum with addition assignment
      //=====================================================================================

      // Conjugate transpose minimum with addition assignment with the given vectors
      {
         test_  = "Conjugate transpose minimum with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initTransposeResults();
            tdres_   += ctrans( min( lhs_, rhs_ ) );
            tsres_   += ctrans( min( lhs_, rhs_ ) );
            trefres_ += ctrans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   += ctrans( min( tlhs_, trhs_ ) );
            sres_   += ctrans( min( tlhs_, trhs_ ) );
            refres_ += ctrans( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Conjugate transpose minimum with addition assignment with evaluated vectors
      {
         test_  = "Conjugate transpose minimum with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initTransposeResults();
            tdres_   += ctrans( min( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   += ctrans( min( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ += ctrans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   += ctrans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   += ctrans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ += ctrans( eval( tref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Conjugate transpose minimum with subtraction assignment
      //=====================================================================================

      // Conjugate transpose minimum with subtraction assignment with the given vectors
      {
         test_  = "Conjugate transpose minimum with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initTransposeResults();
            tdres_   -= ctrans( min( lhs_, rhs_ ) );
            tsres_   -= ctrans( min( lhs_, rhs_ ) );
            trefres_ -= ctrans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   -= ctrans( min( tlhs_, trhs_ ) );
            sres_   -= ctrans( min( tlhs_, trhs_ ) );
            refres_ -= ctrans( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Conjugate transpose minimum with subtraction assignment with evaluated vectors
      {
         test_  = "Conjugate transpose minimum with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initTransposeResults();
            tdres_   -= ctrans( min( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   -= ctrans( min( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ -= ctrans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   -= ctrans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   -= ctrans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ -= ctrans( eval( tref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Conjugate transpose minimum with multiplication assignment
      //=====================================================================================

      // Conjugate transpose minimum with multiplication assignment with the given vectors
      {
         test_  = "Conjugate transpose minimum with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initTransposeResults();
            tdres_   *= ctrans( min( lhs_, rhs_ ) );
            tsres_   *= ctrans( min( lhs_, rhs_ ) );
            trefres_ *= ctrans( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   *= ctrans( min( tlhs_, trhs_ ) );
            sres_   *= ctrans( min( tlhs_, trhs_ ) );
            refres_ *= ctrans( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }

      // Conjugate transpose minimum with multiplication assignment with evaluated vectors
      {
         test_  = "Conjugate transpose minimum with multiplication assignment with evaluated vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initTransposeResults();
            tdres_   *= ctrans( min( eval( lhs_ ), eval( rhs_ ) ) );
            tsres_   *= ctrans( min( eval( lhs_ ), eval( rhs_ ) ) );
            trefres_ *= ctrans( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkTransposeResults<VT1,VT2>();

         try {
            initResults();
            dres_   *= ctrans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            sres_   *= ctrans( min( eval( tlhs_ ), eval( trhs_ ) ) );
            refres_ *= ctrans( eval( tref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Conjugate transpose minimum with division assignment
      //=====================================================================================

      if( blaze::isDivisor( min( lhs_, rhs_ ) ) )
      {
         // Conjugate transpose minimum with division assignment with the given vectors
         {
            test_  = "Conjugate transpose minimum with division assignment with the given vectors";
            error_ = "Failed division assignment operation";

            try {
               initTransposeResults();
               tdres_   /= ctrans( min( lhs_, rhs_ ) );
               tsres_   /= ctrans( min( lhs_, rhs_ ) );
               trefres_ /= ctrans( ref_ );
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkTransposeResults<VT1,VT2>();

            try {
               initResults();
               dres_   /= ctrans( min( tlhs_, trhs_ ) );
               sres_   /= ctrans( min( tlhs_, trhs_ ) );
               refres_ /= ctrans( tref_ );
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkResults<TVT1,TVT2>();
         }

         // Conjugate transpose minimum with division assignment with evaluated vectors
         {
            test_  = "Conjugate transpose minimum with division assignment with evaluated vectors";
            error_ = "Failed division assignment operation";

            try {
               initTransposeResults();
               tdres_   /= ctrans( min( eval( lhs_ ), eval( rhs_ ) ) );
               tsres_   /= ctrans( min( eval( lhs_ ), eval( rhs_ ) ) );
               trefres_ /= ctrans( eval( ref_ ) );
            }
            catch( std::exception& ex ) {
               convertException<VT1,VT2>( ex );
            }

            checkTransposeResults<VT1,VT2>();

            try {
               initResults();
               dres_   /= ctrans( min( eval( tlhs_ ), eval( trhs_ ) ) );
               sres_   /= ctrans( min( eval( tlhs_ ), eval( trhs_ ) ) );
               refres_ /= ctrans( eval( tref_ ) );
            }
            catch( std::exception& ex ) {
               convertException<TVT1,TVT2>( ex );
            }

            checkResults<TVT1,TVT2>();
         }
      }
   }
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the abs dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the abs vector minimum with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case any
// error resulting from the minimum operation or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
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
/*!\brief Testing the conjugate dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the conjugate vector minimum with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case any
// error resulting from the minimum operation or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
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
/*!\brief Testing the \a real dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the \a real vector minimum with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case any
// error resulting from the minimum operation or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
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
/*!\brief Testing the \a imag dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the \a imag vector minimum with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case any
// error resulting from the minimum operation or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
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
/*!\brief Testing the evaluated dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the evaluated vector minimum with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case any
// error resulting from the minimum operation or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
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
/*!\brief Testing the serialized dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the serialized vector minimum with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case any
// error resulting from the minimum operation or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
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
/*!\brief Testing the non-aliased dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the non-aliased vector minimum with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case any
// error resulting from the minimum operation or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
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
/*!\brief Testing the non-SIMD dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the non-SIMD vector minimum with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment. In case any
// error resulting from the minimum operation or the subsequent assignment is detected, a
// \a std::runtime_error exception is thrown.
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
/*!\brief Testing the subvector-wise dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the subvector-wise vector minimum with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment.
// In case any error resulting from the minimum operation or the subsequent assignment is
// detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testSubvectorOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION
   if( BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION > 1 )
   {
      if( lhs_.size() == 0UL )
         return;


      //=====================================================================================
      // Subvector-wise minimum
      //=====================================================================================

      // Subvector-wise minimum with the given vectors
      {
         test_  = "Subvector-wise minimum with the given vectors";
         error_ = "Failed minimum operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size() - index );
               subvector( dres_  , index, size ) = subvector( min( lhs_, rhs_ ), index, size );
               subvector( sres_  , index, size ) = subvector( min( lhs_, rhs_ ), index, size );
               subvector( refres_, index, size ) = subvector( ref_             , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size() - index );
               subvector( tdres_  , index, size ) = subvector( min( tlhs_, trhs_ ), index, size );
               subvector( tsres_  , index, size ) = subvector( min( tlhs_, trhs_ ), index, size );
               subvector( trefres_, index, size ) = subvector( tref_              , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subvector-wise minimum with evaluated vectors
      {
         test_  = "Subvector-wise minimum with evaluated vectors";
         error_ = "Failed minimum operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size() - index );
               subvector( dres_  , index, size ) = subvector( min( eval( lhs_ ), eval( rhs_ ) ), index, size );
               subvector( sres_  , index, size ) = subvector( min( eval( lhs_ ), eval( rhs_ ) ), index, size );
               subvector( refres_, index, size ) = subvector( eval( ref_ )                     , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size() - index );
               subvector( tdres_  , index, size ) = subvector( min( eval( tlhs_ ), eval( trhs_ ) ), index, size );
               subvector( tsres_  , index, size ) = subvector( min( eval( tlhs_ ), eval( trhs_ ) ), index, size );
               subvector( trefres_, index, size ) = subvector( eval( tref_ )                      , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Subvector-wise minimum with addition assignment
      //=====================================================================================

      // Subvector-wise minimum with addition assignment with the given vectors
      {
         test_  = "Subvector-wise minimum with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size() - index );
               subvector( dres_  , index, size ) += subvector( min( lhs_, rhs_ ), index, size );
               subvector( sres_  , index, size ) += subvector( min( lhs_, rhs_ ), index, size );
               subvector( refres_, index, size ) += subvector( ref_             , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size() - index );
               subvector( tdres_  , index, size ) += subvector( min( tlhs_, trhs_ ), index, size );
               subvector( tsres_  , index, size ) += subvector( min( tlhs_, trhs_ ), index, size );
               subvector( trefres_, index, size ) += subvector( tref_              , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subvector-wise minimum with addition assignment with evaluated vectors
      {
         test_  = "Subvector-wise minimum with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size() - index );
               subvector( dres_  , index, size ) += subvector( min( eval( lhs_ ), eval( rhs_ ) ), index, size );
               subvector( sres_  , index, size ) += subvector( min( eval( lhs_ ), eval( rhs_ ) ), index, size );
               subvector( refres_, index, size ) += subvector( eval( ref_ )                     , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size() - index );
               subvector( tdres_  , index, size ) += subvector( min( eval( tlhs_ ), eval( trhs_ ) ), index, size );
               subvector( tsres_  , index, size ) += subvector( min( eval( tlhs_ ), eval( trhs_ ) ), index, size );
               subvector( trefres_, index, size ) += subvector( eval( tref_ )                      , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Subvector-wise minimum with subtraction assignment
      //=====================================================================================

      // Subvector-wise minimum with subtraction assignment with the given vectors
      {
         test_  = "Subvector-wise minimum with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size() - index );
               subvector( dres_  , index, size ) -= subvector( min( lhs_, rhs_ ), index, size );
               subvector( sres_  , index, size ) -= subvector( min( lhs_, rhs_ ), index, size );
               subvector( refres_, index, size ) -= subvector( ref_             , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size() - index );
               subvector( tdres_  , index, size ) -= subvector( min( tlhs_, trhs_ ), index, size );
               subvector( tsres_  , index, size ) -= subvector( min( tlhs_, trhs_ ), index, size );
               subvector( trefres_, index, size ) -= subvector( tref_              , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subvector-wise minimum with subtraction assignment with evaluated vectors
      {
         test_  = "Subvector-wise minimum with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size() - index );
               subvector( dres_  , index, size ) -= subvector( min( eval( lhs_ ), eval( rhs_ ) ), index, size );
               subvector( sres_  , index, size ) -= subvector( min( eval( lhs_ ), eval( rhs_ ) ), index, size );
               subvector( refres_, index, size ) -= subvector( eval( ref_ )                     , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size() - index );
               subvector( tdres_  , index, size ) -= subvector( min( eval( tlhs_ ), eval( trhs_ ) ), index, size );
               subvector( tsres_  , index, size ) -= subvector( min( eval( tlhs_ ), eval( trhs_ ) ), index, size );
               subvector( trefres_, index, size ) -= subvector( eval( tref_ )                      , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Subvector-wise minimum with multiplication assignment
      //=====================================================================================

      // Subvector-wise minimum with multiplication assignment with the given vectors
      {
         test_  = "Subvector-wise minimum with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size() - index );
               subvector( dres_  , index, size ) *= subvector( min( lhs_, rhs_ ), index, size );
               subvector( sres_  , index, size ) *= subvector( min( lhs_, rhs_ ), index, size );
               subvector( refres_, index, size ) *= subvector( ref_             , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size() - index );
               subvector( tdres_  , index, size ) *= subvector( min( tlhs_, trhs_ ), index, size );
               subvector( tsres_  , index, size ) *= subvector( min( tlhs_, trhs_ ), index, size );
               subvector( trefres_, index, size ) *= subvector( tref_              , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subvector-wise minimum with multiplication assignment with evaluated vectors
      {
         test_  = "Subvector-wise minimum with multiplication assignment with evaluated vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size() - index );
               subvector( dres_  , index, size ) *= subvector( min( eval( lhs_ ), eval( rhs_ ) ), index, size );
               subvector( sres_  , index, size ) *= subvector( min( eval( lhs_ ), eval( rhs_ ) ), index, size );
               subvector( refres_, index, size ) *= subvector( eval( ref_ )                     , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size() - index );
               subvector( tdres_  , index, size ) *= subvector( min( eval( tlhs_ ), eval( trhs_ ) ), index, size );
               subvector( tsres_  , index, size ) *= subvector( min( eval( tlhs_ ), eval( trhs_ ) ), index, size );
               subvector( trefres_, index, size ) *= subvector( eval( tref_ )                      , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Subvector-wise minimum with division assignment
      //=====================================================================================

      // Subvector-wise minimum with division assignment with the given vectors
      {
         test_  = "Subvector-wise minimum with division assignment with the given vectors";
         error_ = "Failed division assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size() - index );
               if( !blaze::isDivisor( subvector( min( lhs_, rhs_ ), index, size ) ) ) continue;
               subvector( dres_  , index, size ) /= subvector( min( lhs_, rhs_ ), index, size );
               subvector( sres_  , index, size ) /= subvector( min( lhs_, rhs_ ), index, size );
               subvector( refres_, index, size ) /= subvector( ref_             , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size() - index );
               if( !blaze::isDivisor( subvector( min( tlhs_, trhs_ ), index, size ) ) ) continue;
               subvector( tdres_  , index, size ) /= subvector( min( tlhs_, trhs_ ), index, size );
               subvector( tsres_  , index, size ) /= subvector( min( tlhs_, trhs_ ), index, size );
               subvector( trefres_, index, size ) /= subvector( tref_              , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Subvector-wise minimum with division assignment with evaluated vectors
      {
         test_  = "Subvector-wise minimum with division assignment with evaluated vectors";
         error_ = "Failed division assignment operation";

         try {
            initResults();
            for( size_t index=0UL, size=0UL; index<lhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, lhs_.size() - index );
               if( !blaze::isDivisor( subvector( min( lhs_, rhs_ ), index, size ) ) ) continue;
               subvector( dres_  , index, size ) /= subvector( min( eval( lhs_ ), eval( rhs_ ) ), index, size );
               subvector( sres_  , index, size ) /= subvector( min( eval( lhs_ ), eval( rhs_ ) ), index, size );
               subvector( refres_, index, size ) /= subvector( eval( ref_ )                     , index, size );
            }
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            for( size_t index=0UL, size=0UL; index<tlhs_.size(); index+=size ) {
               size = blaze::rand<size_t>( 1UL, tlhs_.size() - index );
               if( !blaze::isDivisor( subvector( min( tlhs_, trhs_ ), index, size ) ) ) continue;
               subvector( tdres_  , index, size ) /= subvector( min( eval( tlhs_ ), eval( trhs_ ) ), index, size );
               subvector( tsres_  , index, size ) /= subvector( min( eval( tlhs_ ), eval( trhs_ ) ), index, size );
               subvector( trefres_, index, size ) /= subvector( eval( tref_ )                      , index, size );
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
         auto sv = subvector( min( lhs_, rhs_ ), 1UL, lhs_.size() );

         std::ostringstream oss;
         oss << " Test: Subvector construction\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
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
         auto sv = subvector( min( lhs_, rhs_ ), lhs_.size(), 1UL );

         std::ostringstream oss;
         oss << " Test: Subvector construction\n"
             << " Error: Setup of out-of-bounds subvector succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
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
/*!\brief Skipping the subvector-wise dense vector/dense vector minimum operation.
//
// \return void
//
// This function is called in case the subvector-wise vector/vector minimum operation is not
// available for the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testSubvectorOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the elements-wise dense vector/dense vector minimum operation.
//
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the elements-wise vector minimum with plain assignment, addition
// assignment, subtraction assignment, multiplication assignment, and division assignment.
// In case any error resulting from the minimum operation or the subsequent assignment is
// detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testElementsOperation( blaze::TrueType )
{
#if BLAZETEST_MATHTEST_TEST_ELEMENTS_OPERATION
   if( BLAZETEST_MATHTEST_TEST_ELEMENTS_OPERATION > 1 )
   {
      if( lhs_.size() == 0UL )
         return;


      std::vector<size_t> indices( lhs_.size() );
      std::iota( indices.begin(), indices.end(), 0UL );
      std::random_shuffle( indices.begin(), indices.end() );


      //=====================================================================================
      // Elements-wise minimum
      //=====================================================================================

      // Elements-wise minimum with the given vectors
      {
         test_  = "Elements-wise minimum with the given vectors";
         error_ = "Failed minimum operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) = elements( min( lhs_, rhs_ ), &indices[index], n );
               elements( sres_  , &indices[index], n ) = elements( min( lhs_, rhs_ ), &indices[index], n );
               elements( refres_, &indices[index], n ) = elements( ref_             , &indices[index], n );
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
               elements( tdres_  , &indices[index], n ) = elements( min( tlhs_, trhs_ ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) = elements( min( tlhs_, trhs_ ), &indices[index], n );
               elements( trefres_, &indices[index], n ) = elements( tref_              , &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Elements-wise minimum with evaluated vectors
      {
         test_  = "Elements-wise minimum with evaluated vectors";
         error_ = "Failed minimum operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) = elements( min( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               elements( sres_  , &indices[index], n ) = elements( min( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               elements( refres_, &indices[index], n ) = elements( eval( ref_ )                     , &indices[index], n );
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
               elements( tdres_  , &indices[index], n ) = elements( min( eval( tlhs_ ), eval( trhs_ ) ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) = elements( min( eval( tlhs_ ), eval( trhs_ ) ), &indices[index], n );
               elements( trefres_, &indices[index], n ) = elements( eval( tref_ )                      , &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Elements-wise minimum with addition assignment
      //=====================================================================================

      // Elements-wise minimum with addition assignment with the given vectors
      {
         test_  = "Elements-wise minimum with addition assignment with the given vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) += elements( min( lhs_, rhs_ ), &indices[index], n );
               elements( sres_  , &indices[index], n ) += elements( min( lhs_, rhs_ ), &indices[index], n );
               elements( refres_, &indices[index], n ) += elements( ref_             , &indices[index], n );
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
               elements( tdres_  , &indices[index], n ) += elements( min( tlhs_, trhs_ ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) += elements( min( tlhs_, trhs_ ), &indices[index], n );
               elements( trefres_, &indices[index], n ) += elements( tref_              , &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Elements-wise minimum with addition assignment with evaluated vectors
      {
         test_  = "Elements-wise minimum with addition assignment with evaluated vectors";
         error_ = "Failed addition assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) += elements( min( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               elements( sres_  , &indices[index], n ) += elements( min( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               elements( refres_, &indices[index], n ) += elements( eval( ref_ )                     , &indices[index], n );
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
               elements( tdres_  , &indices[index], n ) += elements( min( eval( tlhs_ ), eval( trhs_ ) ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) += elements( min( eval( tlhs_ ), eval( trhs_ ) ), &indices[index], n );
               elements( trefres_, &indices[index], n ) += elements( eval( tref_ )                      , &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Elements-wise minimum with subtraction assignment
      //=====================================================================================

      // Elements-wise minimum with subtraction assignment with the given vectors
      {
         test_  = "Elements-wise minimum with subtraction assignment with the given vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) -= elements( min( lhs_, rhs_ ), &indices[index], n );
               elements( sres_  , &indices[index], n ) -= elements( min( lhs_, rhs_ ), &indices[index], n );
               elements( refres_, &indices[index], n ) -= elements( ref_             , &indices[index], n );
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
               elements( tdres_  , &indices[index], n ) -= elements( min( tlhs_, trhs_ ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) -= elements( min( tlhs_, trhs_ ), &indices[index], n );
               elements( trefres_, &indices[index], n ) -= elements( tref_              , &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Elements-wise minimum with subtraction assignment with evaluated vectors
      {
         test_  = "Elements-wise minimum with subtraction assignment with evaluated vectors";
         error_ = "Failed subtraction assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) -= elements( min( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               elements( sres_  , &indices[index], n ) -= elements( min( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               elements( refres_, &indices[index], n ) -= elements( eval( ref_ )                     , &indices[index], n );
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
               elements( tdres_  , &indices[index], n ) -= elements( min( eval( tlhs_ ), eval( trhs_ ) ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) -= elements( min( eval( tlhs_ ), eval( trhs_ ) ), &indices[index], n );
               elements( trefres_, &indices[index], n ) -= elements( eval( tref_ )                      , &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Elements-wise minimum with multiplication assignment
      //=====================================================================================

      // Elements-wise minimum with multiplication assignment with the given vectors
      {
         test_  = "Elements-wise minimum with multiplication assignment with the given vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) *= elements( min( lhs_, rhs_ ), &indices[index], n );
               elements( sres_  , &indices[index], n ) *= elements( min( lhs_, rhs_ ), &indices[index], n );
               elements( refres_, &indices[index], n ) *= elements( ref_             , &indices[index], n );
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
               elements( tdres_  , &indices[index], n ) *= elements( min( tlhs_, trhs_ ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) *= elements( min( tlhs_, trhs_ ), &indices[index], n );
               elements( trefres_, &indices[index], n ) *= elements( tref_              , &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Elements-wise minimum with multiplication assignment with evaluated vectors
      {
         test_  = "Elements-wise minimum with multiplication assignment with evaluated vectors";
         error_ = "Failed multiplication assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               elements( dres_  , &indices[index], n ) *= elements( min( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               elements( sres_  , &indices[index], n ) *= elements( min( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               elements( refres_, &indices[index], n ) *= elements( eval( ref_ )                     , &indices[index], n );
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
               elements( tdres_  , &indices[index], n ) *= elements( min( eval( tlhs_ ), eval( trhs_ ) ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) *= elements( min( eval( tlhs_ ), eval( trhs_ ) ), &indices[index], n );
               elements( trefres_, &indices[index], n ) *= elements( eval( tref_ )                      , &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }


      //=====================================================================================
      // Elements-wise minimum with division assignment
      //=====================================================================================

      // Elements-wise minimum with division assignment with the given vectors
      {
         test_  = "Elements-wise minimum with division assignment with the given vectors";
         error_ = "Failed division assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               if( !blaze::isDivisor( elements( min( lhs_, rhs_ ), &indices[index], n ) ) ) continue;
               elements( dres_  , &indices[index], n ) /= elements( min( lhs_, rhs_ ), &indices[index], n );
               elements( sres_  , &indices[index], n ) /= elements( min( lhs_, rhs_ ), &indices[index], n );
               elements( refres_, &indices[index], n ) /= elements( ref_             , &indices[index], n );
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
               if( !blaze::isDivisor( elements( min( tlhs_, trhs_ ), &indices[index], n ) ) ) continue;
               elements( tdres_  , &indices[index], n ) /= elements( min( tlhs_, trhs_ ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) /= elements( min( tlhs_, trhs_ ), &indices[index], n );
               elements( trefres_, &indices[index], n ) /= elements( tref_              , &indices[index], n );
            }
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Elements-wise minimum with division assignment with evaluated vectors
      {
         test_  = "Elements-wise minimum with division assignment with evaluated vectors";
         error_ = "Failed division assignment operation";

         try {
            initResults();
            for( size_t index=0UL, n=0UL; index<indices.size(); index+=n ) {
               n = blaze::rand<size_t>( 1UL, indices.size() - index );
               if( !blaze::isDivisor( elements( min( lhs_, rhs_ ), &indices[index], n ) ) ) continue;
               elements( dres_  , &indices[index], n ) /= elements( min( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               elements( sres_  , &indices[index], n ) /= elements( min( eval( lhs_ ), eval( rhs_ ) ), &indices[index], n );
               elements( refres_, &indices[index], n ) /= elements( eval( ref_ )                     , &indices[index], n );
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
               if( !blaze::isDivisor( elements( min( tlhs_, trhs_ ), &indices[index], n ) ) ) continue;
               elements( tdres_  , &indices[index], n ) /= elements( min( eval( tlhs_ ), eval( trhs_ ) ), &indices[index], n );
               elements( tsres_  , &indices[index], n ) /= elements( min( eval( tlhs_ ), eval( trhs_ ) ), &indices[index], n );
               elements( trefres_, &indices[index], n ) /= elements( eval( tref_ )                      , &indices[index], n );
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
         auto e = elements( min( lhs_, rhs_ ), blaze::index_sequence<128>() );

         std::ostringstream oss;
         oss << " Test: Elements construction (index_sequence)\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
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
         auto e = elements( min( lhs_, rhs_ ), { lhs_.size() } );

         std::ostringstream oss;
         oss << " Test: Elements construction (initializer_list)\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
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
         auto e = elements( min( lhs_, rhs_ ), [index=lhs_.size()]( size_t ){ return index; }, 1UL );

         std::ostringstream oss;
         oss << " Test: Elements construction (lambda)\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Random seed = " << blaze::getSeed() << "\n"
             << "   Left-hand side dense vector type:\n"
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
/*!\brief Skipping the elements-wise dense vector/dense vector minimum operation.
//
// \return void
//
// This function is called in case the elements-wise vector/vector minimum operation is not
// available for the given vector types \a VT1 and \a VT2.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::testElementsOperation( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the customized dense vector/dense vector minimum operation.
//
// \param op The custom operation to be tested.
// \param name The human-readable name of the operation.
// \return void
// \exception std::runtime_error Minimum error detected.
//
// This function tests the vector minimum with plain assignment, addition assignment,
// subtraction assignment, multiplication assignment, and division assignment in combination with
// a custom operation. In case any error resulting from the minimum operation or the subsequent
// assignment is detected, a \a std::runtime_error exception is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
template< typename OP >   // Type of the custom operation
void OperationTest<VT1,VT2>::testCustomOperation( OP op, const std::string& name )
{
   //=====================================================================================
   // Customized minimum
   //=====================================================================================

   // Customized minimum with the given vectors
   {
      test_  = "Customized minimum with the given vectors (" + name + ")";
      error_ = "Failed minimum operation";

      try {
         initResults();
         dres_   = op( min( lhs_, rhs_ ) );
         sres_   = op( min( lhs_, rhs_ ) );
         refres_ = op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   = op( min( tlhs_, trhs_ ) );
         tsres_   = op( min( tlhs_, trhs_ ) );
         trefres_ = op( tref_ );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }

   // Customized minimum with evaluated vectors
   {
      test_  = "Customized minimum with evaluated vectors (" + name + ")";
      error_ = "Failed minimum operation";

      try {
         initResults();
         dres_   = op( min( eval( lhs_ ), eval( rhs_ ) ) );
         sres_   = op( min( eval( lhs_ ), eval( rhs_ ) ) );
         refres_ = op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   = op( min( eval( tlhs_ ), eval( trhs_ ) ) );
         tsres_   = op( min( eval( tlhs_ ), eval( trhs_ ) ) );
         trefres_ = op( eval( tref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }


   //=====================================================================================
   // Customized minimum with addition assignment
   //=====================================================================================

   // Customized minimum with addition assignment with the given vectors
   {
      test_  = "Customized minimum with addition assignment with the given vectors (" + name + ")";
      error_ = "Failed addition assignment operation";

      try {
         initResults();
         dres_   += op( min( lhs_, rhs_ ) );
         sres_   += op( min( lhs_, rhs_ ) );
         refres_ += op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   += op( min( tlhs_, trhs_ ) );
         tsres_   += op( min( tlhs_, trhs_ ) );
         trefres_ += op( tref_ );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }

   // Customized minimum with addition assignment with evaluated vectors
   {
      test_  = "Customized minimum with addition assignment with evaluated vectors (" + name + ")";
      error_ = "Failed addition assignment operation";

      try {
         initResults();
         dres_   += op( min( eval( lhs_ ), eval( rhs_ ) ) );
         sres_   += op( min( eval( lhs_ ), eval( rhs_ ) ) );
         refres_ += op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   += op( min( eval( tlhs_ ), eval( trhs_ ) ) );
         tsres_   += op( min( eval( tlhs_ ), eval( trhs_ ) ) );
         trefres_ += op( eval( tref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }


   //=====================================================================================
   // Customized minimum with subtraction assignment
   //=====================================================================================

   // Customized minimum with subtraction assignment with the given vectors
   {
      test_  = "Customized minimum with subtraction assignment with the given vectors (" + name + ")";
      error_ = "Failed subtraction assignment operation";

      try {
         initResults();
         dres_   -= op( min( lhs_, rhs_ ) );
         sres_   -= op( min( lhs_, rhs_ ) );
         refres_ -= op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   -= op( min( tlhs_, trhs_ ) );
         tsres_   -= op( min( tlhs_, trhs_ ) );
         trefres_ -= op( tref_ );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }

   // Customized minimum with subtraction assignment with evaluated vectors
   {
      test_  = "Customized minimum with subtraction assignment with evaluated vectors (" + name + ")";
      error_ = "Failed subtraction assignment operation";

      try {
         initResults();
         dres_   -= op( min( eval( lhs_ ), eval( rhs_ ) ) );
         sres_   -= op( min( eval( lhs_ ), eval( rhs_ ) ) );
         refres_ -= op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   -= op( min( eval( tlhs_ ), eval( trhs_ ) ) );
         tsres_   -= op( min( eval( tlhs_ ), eval( trhs_ ) ) );
         trefres_ -= op( eval( tref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }


   //=====================================================================================
   // Customized minimum with multiplication assignment
   //=====================================================================================

   // Customized minimum with multiplication assignment with the given vectors
   {
      test_  = "Customized minimum with multiplication assignment with the given vectors (" + name + ")";
      error_ = "Failed multiplication assignment operation";

      try {
         initResults();
         dres_   *= op( min( lhs_, rhs_ ) );
         sres_   *= op( min( lhs_, rhs_ ) );
         refres_ *= op( ref_ );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   *= op( min( tlhs_, trhs_ ) );
         tsres_   *= op( min( tlhs_, trhs_ ) );
         trefres_ *= op( tref_ );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }

   // Customized minimum with multiplication assignment with evaluated vectors
   {
      test_  = "Customized minimum with multiplication assignment with evaluated vectors (" + name + ")";
      error_ = "Failed multiplication assignment operation";

      try {
         initResults();
         dres_   *= op( min( eval( lhs_ ), eval( rhs_ ) ) );
         sres_   *= op( min( eval( lhs_ ), eval( rhs_ ) ) );
         refres_ *= op( eval( ref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<VT1,VT2>( ex );
      }

      checkResults<VT1,VT2>();

      try {
         initTransposeResults();
         tdres_   *= op( min( eval( tlhs_ ), eval( trhs_ ) ) );
         tsres_   *= op( min( eval( tlhs_ ), eval( trhs_ ) ) );
         trefres_ *= op( eval( tref_ ) );
      }
      catch( std::exception& ex ) {
         convertException<TVT1,TVT2>( ex );
      }

      checkTransposeResults<TVT1,TVT2>();
   }


   //=====================================================================================
   // Customized minimum with division assignment
   //=====================================================================================

   if( blaze::isDivisor( op( min( lhs_, rhs_ ) ) ) )
   {
      // Customized minimum with division assignment with the given vectors
      {
         test_  = "Customized minimum with division assignment with the given vectors (" + name + ")";
         error_ = "Failed division assignment operation";

         try {
            initResults();
            dres_   /= op( min( lhs_, rhs_ ) );
            sres_   /= op( min( lhs_, rhs_ ) );
            refres_ /= op( ref_ );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   /= op( min( tlhs_, trhs_ ) );
            tsres_   /= op( min( tlhs_, trhs_ ) );
            trefres_ /= op( tref_ );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
      }

      // Customized minimum with division assignment with evaluated vectors
      {
         test_  = "Customized minimum with division assignment with evaluated vectors (" + name + ")";
         error_ = "Failed division assignment operation";

         try {
            initResults();
            dres_   /= op( min( eval( lhs_ ), eval( rhs_ ) ) );
            sres_   /= op( min( eval( lhs_ ), eval( rhs_ ) ) );
            refres_ /= op( eval( ref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<VT1,VT2>( ex );
         }

         checkResults<VT1,VT2>();

         try {
            initTransposeResults();
            tdres_   /= op( min( eval( tlhs_ ), eval( trhs_ ) ) );
            tsres_   /= op( min( eval( tlhs_ ), eval( trhs_ ) ) );
            trefres_ /= op( eval( tref_ ) );
         }
         catch( std::exception& ex ) {
            convertException<TVT1,TVT2>( ex );
         }

         checkTransposeResults<TVT1,TVT2>();
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
// This function is called after each test case to check and compare the computed results. The
// two template arguments \a LT and \a RT indicate the types of the left-hand side and right-hand
// side operands used for the computations.
*/
template< typename VT1    // Type of the left-hand side dense vector
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
          << "   Left-hand side dense " << ( IsRowVector<LT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
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
          << "   Left-hand side dense " << ( IsRowVector<LT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
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
template< typename VT1    // Type of the left-hand side dense vector
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
          << "   Left-hand side dense " << ( IsRowVector<LT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
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
          << "   Left-hand side dense " << ( IsRowVector<LT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
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
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::initResults()
{
   const blaze::UnderlyingBuiltin_t<DRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<DRE> max( randmax );

   resize( dres_, size( lhs_ ) );
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
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void OperationTest<VT1,VT2>::initTransposeResults()
{
   const blaze::UnderlyingBuiltin_t<TDRE> min( randmin );
   const blaze::UnderlyingBuiltin_t<TDRE> max( randmax );

   resize( tdres_, size( tlhs_ ) );
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
template< typename VT1    // Type of the left-hand side dense vector
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
       << "   Left-hand side dense " << ( IsRowVector<LT>::value ? ( "row" ) : ( "column" ) ) << " vector type:\n"
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
/*!\brief Testing the vector minimum operation between two specific vector types.
//
// \param creator1 The creator for the left-hand side dense vector.
// \param creator2 The creator for the right-hand side dense vector.
// \return void
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
void runTest( const Creator<VT1>& creator1, const Creator<VT2>& creator2 )
{
#if BLAZETEST_MATHTEST_TEST_MINIMUM
   if( BLAZETEST_MATHTEST_TEST_MINIMUM > 1 )
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
/*!\brief Macro for the definition of a dense vector/dense vector minimum test case.
*/
#define DEFINE_DVECDVECMIN_OPERATION_TEST( VT1, VT2 ) \
   extern template class blazetest::mathtest::operations::dvecdvecmin::OperationTest<VT1,VT2>
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of a dense vector/dense vector minimum test case.
*/
#define RUN_DVECDVECMIN_OPERATION_TEST( C1, C2 ) \
   blazetest::mathtest::operations::dvecdvecmin::runTest( C1, C2 )
/*! \endcond */
//*************************************************************************************************

} // namespace dvecdvecmin

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
