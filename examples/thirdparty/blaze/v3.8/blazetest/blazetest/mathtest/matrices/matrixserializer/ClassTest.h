//=================================================================================================
/*!
//  \file blazetest/mathtest/matrices/matrixserializer/ClassTest.h
//  \brief Header file for the MatrixSerializer class test
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

#ifndef _BLAZETEST_MATHTEST_MATRICES_MATRIXSERIALIZER_CLASSTEST_H_
#define _BLAZETEST_MATHTEST_MATRICES_MATRIXSERIALIZER_CLASSTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/Aliases.h>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/constraints/Matrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/Submatrix.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/util/Random.h>
#include <blaze/util/serialization/Archive.h>


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace matrixserializer {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the MatrixSerializer class.
//
// This class represents a test suite for the blaze::MatrixSerializer class. It performs a
// series of runtime tests with different matrix types to test the serialization of both
// dense and sparse matrices.
*/
class ClassTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit ClassTest();
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
   void testEmptyMatrices ();
   void testRandomMatrices();
   void testFailures      ();

   template< size_t M, size_t N, typename MT >
   void runAllTests( const MT& src );

   template< size_t M, size_t N, typename MT >
   void runStaticMatrixTests( const MT& src );

   template< typename MT >
   void runDynamicMatrixTests( const MT& src );

   template< size_t M, size_t N, typename MT >
   void runDenseSubmatrixTests( const MT& src );

   template< typename MT >
   void runCompressedMatrixTests( const MT& src );

   template< size_t M, size_t N, typename MT >
   void runSparseSubmatrixTests( const MT& src );

   template< typename MT1, typename MT2 >
   void runTest( const MT1& src, MT2& dst );

   template< typename Archive, typename MT >
   void testSerialization( Archive& archive, const MT& src );

   template< typename Archive, typename MT >
   void testDeserialization( Archive& archive, MT& dst );

   template< typename MT1, typename MT2 >
   void compareMatrices( const MT1& src, const MT2& dst );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Execution of several (de-)serialization tests with the given source matrix.
//
// \param src The source matrix to be tested.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the matrix (de-)serialization with the given matrix. The matrix is
// serialized and deserialized several times, using instances of StaticMatrix, DynamicMatrix,
// and CompressedMatrix as destination matrix type. In case an error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< size_t M       // Number of rows of the matrix
        , size_t N       // Number of columns of the matrix
        , typename MT >  // Type of the matrix
void ClassTest::runAllTests( const MT& src )
{
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( MT );

   runStaticMatrixTests<M,N>   ( src );
   runDynamicMatrixTests       ( src );
   runDenseSubmatrixTests<M,N> ( src );
   runCompressedMatrixTests    ( src );
   runSparseSubmatrixTests<M,N>( src );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Execution of several (de-)serialization tests with the given source matrix.
//
// \param src The source matrix to be tested.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the matrix (de-)serialization with the given matrix. The matrix is
// serialized and deserialized several times, using instances of StaticMatrix as destination
// matrix type. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t M       // Number of rows of the matrix
        , size_t N       // Number of columns of the matrix
        , typename MT >  // Type of the matrix
void ClassTest::runStaticMatrixTests( const MT& src )
{
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( MT );

   using ET = blaze::ElementType_t<MT>;

   {
      blaze::StaticMatrix<ET,M,N,blaze::rowMajor> dst;
      randomize( dst );
      runTest( src, dst );
   }

   {
      blaze::StaticMatrix<ET,M,N,blaze::columnMajor> dst;
      randomize( dst );
      runTest( src, dst );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Execution of several (de-)serialization tests with the given source matrix.
//
// \param src The source matrix to be tested.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the matrix (de-)serialization with the given matrix. The matrix is
// serialized and deserialized several times, using instances of DynamicMatrix as destination
// matrix type. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT >  // Type of the matrix
void ClassTest::runDynamicMatrixTests( const MT& src )
{
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( MT );

   using ET = blaze::ElementType_t<MT>;

   {
      blaze::DynamicMatrix<ET,blaze::rowMajor> dst;
      runTest( src, dst );
   }

   {
      blaze::DynamicMatrix<ET,blaze::columnMajor> dst;
      runTest( src, dst );
   }

   {
      blaze::DynamicMatrix<ET,blaze::rowMajor> dst( 43UL, 37UL );
      randomize( dst );
      runTest( src, dst );
   }

   {
      blaze::DynamicMatrix<ET,blaze::columnMajor> dst( 37UL, 43UL );
      randomize( dst );
      runTest( src, dst );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Execution of several (de-)serialization tests with the given source matrix.
//
// \param src The source matrix to be tested.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the matrix (de-)serialization with the given matrix. The matrix is
// serialized and deserialized several times, using instances of Submatrix as destination
// matrix type. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t M       // Number of rows of the matrix
        , size_t N       // Number of columns of the matrix
        , typename MT >  // Type of the matrix
void ClassTest::runDenseSubmatrixTests( const MT& src )
{
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( MT );

   using ET = blaze::ElementType_t<MT>;
   using RM = blaze::DynamicMatrix<ET,blaze::rowMajor>;
   using CM = blaze::DynamicMatrix<ET,blaze::columnMajor>;

   {
      RM mat( M, N );
      blaze::Submatrix<RM> dst( mat, 0UL, 0UL, M, N );
      randomize( dst );
      runTest( src, dst );
   }

   {
      CM mat( M, N );
      blaze::Submatrix<CM> dst( mat, 0UL, 0UL, M, N );
      randomize( dst );
      runTest( src, dst );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Execution of several (de-)serialization tests with the given source matrix.
//
// \param src The source matrix to be tested.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the matrix (de-)serialization with the given matrix. The matrix is
// serialized and deserialized several times, using instances of CompressedMatrix as destination
// matrix type. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT >  // Type of the matrix
void ClassTest::runCompressedMatrixTests( const MT& src )
{
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( MT );

   using ET = blaze::ElementType_t<MT>;

   {
      blaze::CompressedMatrix<ET,blaze::rowMajor> dst;
      runTest( src, dst );
   }

   {
      blaze::CompressedMatrix<ET,blaze::columnMajor> dst;
      runTest( src, dst );
   }

   {
      blaze::CompressedMatrix<ET,blaze::rowMajor> dst( 43UL, 37UL );
      randomize( dst );
      runTest( src, dst );
   }

   {
      blaze::CompressedMatrix<ET,blaze::columnMajor> dst( 37UL, 43UL );
      randomize( dst );
      runTest( src, dst );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Execution of several (de-)serialization tests with the given source matrix.
//
// \param src The source matrix to be tested.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the matrix (de-)serialization with the given matrix. The matrix is
// serialized and deserialized several times, using instances of Submatrix as destination
// matrix type. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t M       // Number of rows of the matrix
        , size_t N       // Number of columns of the matrix
        , typename MT >  // Type of the matrix
void ClassTest::runSparseSubmatrixTests( const MT& src )
{
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( MT );

   using ET = blaze::ElementType_t<MT>;
   using RM = blaze::CompressedMatrix<ET,blaze::rowMajor>;
   using CM = blaze::CompressedMatrix<ET,blaze::columnMajor>;

   {
      RM mat( M, N );
      blaze::Submatrix<RM> dst( mat, 0UL, 0UL, M, N );
      randomize( dst );
      runTest( src, dst );
   }

   {
      CM mat( M, N );
      blaze::Submatrix<CM> dst( mat, 0UL, 0UL, M, N );
      randomize( dst );
      runTest( src, dst );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Running a single (de-)serialization test with the given pair of matrices.
//
// \param src The source matrix to be serialized.
// \param dst The destination matrix to be reconstituted.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the (de-)serialization process with the given pair of matrices. The
// source matrix is serialized and the destination matrix is reconstituted from the resulting
// archive. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the source matrix
        , typename MT2 >  // Type of the destination matrix
void ClassTest::runTest( const MT1& src, MT2& dst )
{
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( MT2 );

   blaze::Archive<std::stringstream> archive;

   testSerialization  ( archive, src );
   testDeserialization( archive, dst );
   compareMatrices    ( src, dst );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the serialization process.
//
// \param archive The archive to be written.
// \param src The source matrix to be serialized.
// \return void
// \exception std::runtime_error Error detected.
//
// This function test the serialization process with the given archive and matrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Archive  // Type of the archive
        , typename MT >     // Type of the matrix
void ClassTest::testSerialization( Archive& archive, const MT& src )
{
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( MT );

   using blaze::IsDenseMatrix;

   try {
      archive << src;
   }
   catch( std::runtime_error& ex ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Serialization failed\n"
          << " Details:\n"
          << "   " << ( IsDenseMatrix<MT>::value ? ( "Dense" ) : ( "Sparse" ) ) << " matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Matrix:\n" << src << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the deserialization process.
//
// \param archive The archive to be read.
// \param dst The source matrix to be reconstituted.
// \return void
// \exception std::runtime_error Error detected.
//
// This function test the deserialization process with the given archive and matrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Archive  // Type of the archive
        , typename MT >     // Type of the matrix
void ClassTest::testDeserialization( Archive& archive, MT& dst )
{
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( MT );

   using blaze::IsDenseMatrix;

   try {
      archive >> dst;
   }
   catch( std::runtime_error& ex )
   {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Deserialization failed\n"
          << " Details:\n"
          << "   " << ( IsDenseMatrix<MT>::value ? ( "Dense" ) : ( "Sparse" ) ) << " matrix type:\n"
          << "     " << typeid( MT ).name() << "\n"
          << "   Matrix:\n" << dst << "\n"
          << "   Error message: " << ex.what() << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Comparison of a source and destination matrix.
//
// \param src The source matrix.
// \param dst The destination matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a comparison between the given source and destination matrix. In
// case the matrix are not equal, a \a std::runtime_error exception is thrown.
*/
template< typename MT1    // Type of the source matrix
        , typename MT2 >  // Type of the destination matrix
void ClassTest::compareMatrices( const MT1& src, const MT2& dst )
{
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE( MT2 );

   using blaze::IsDenseMatrix;

   if( src != dst ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Matrix comparison failed\n"
          << " Details:\n"
          << "   " << ( IsDenseMatrix<MT1>::value ? ( "Dense" ) : ( "Sparse" ) ) << " source matrix type:\n"
          << "     " << typeid( MT1 ).name() << "\n"
          << "   " << ( IsDenseMatrix<MT2>::value ? ( "Dense" ) : ( "Sparse" ) ) << " destination matrix type:\n"
          << "     " << typeid( MT2 ).name() << "\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the functionality of the MatrixSerializer class.
//
// \return void
*/
void runTest()
{
   ClassTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the MatrixSerializer class test.
*/
#define RUN_MATRIXSERIALIZER_CLASS_TEST \
   blazetest::mathtest::matrices::matrixserializer::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace matrixserializer

} // namespace matrices

} // namespace mathtest

} // namespace blazetest

#endif
