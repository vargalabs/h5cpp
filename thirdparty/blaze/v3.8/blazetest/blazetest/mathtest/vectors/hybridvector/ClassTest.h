//=================================================================================================
/*!
//  \file blazetest/mathtest/vectors/hybridvector/ClassTest.h
//  \brief Header file for the HybridVector class test
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

#ifndef _BLAZETEST_MATHTEST_VECTORS_HYBRIDVECTOR_CLASSTEST_H_
#define _BLAZETEST_MATHTEST_VECTORS_HYBRIDVECTOR_CLASSTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <array>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/HybridVector.h>
#include <blaze/math/simd/SIMDTrait.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/AlignedAllocator.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/typetraits/AlignmentOf.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace vectors {

namespace hybridvector {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the HybridVector class template.
//
// This class represents a test suite for the blaze::HybridVector class template. It performs
// a series of both compile time as well as runtime tests.
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
   template< typename Type >
   void testAlignment( const std::string& type );

   void testConstructors();
   void testAssignment  ();
   void testAddAssign   ();
   void testSubAssign   ();
   void testMultAssign  ();
   void testDivAssign   ();
   void testCrossAssign ();
   void testScaling     ();
   void testSubscript   ();
   void testAt          ();
   void testIterator    ();
   void testNonZeros    ();
   void testReset       ();
   void testClear       ();
   void testResize      ();
   void testExtend      ();
   void testSwap        ();
   void testIsDefault   ();

   template< typename Type >
   void checkSize( const Type& vector, size_t expectedSize ) const;

   template< typename Type >
   void checkCapacity( const Type& vector, size_t minCapacity ) const;

   template< typename Type >
   void checkNonZeros( const Type& vector, size_t nonzeros ) const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   using VT  = blaze::HybridVector<int,4UL,blaze::rowVector>;     //!< Type of the hybrid vector.
   using TVT = blaze::HybridVector<int,4UL,blaze::columnVector>;  //!< Transpose hybrid vector type.

   using RVT  = VT::Rebind<double>::Other;   //!< Rebound hybrid vector type.
   using TRVT = TVT::Rebind<double>::Other;  //!< Transpose rebound hybrid vector type.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT                  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT::ResultType      );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT::TransposeType   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( TVT                 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( TVT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( TVT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RVT                 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RVT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RVT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( TRVT                );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( TRVT::ResultType    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( TRVT::TransposeType );

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( VT                  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( VT::ResultType      );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( VT::TransposeType   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( TVT                 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( TVT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( TVT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( RVT                 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( RVT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( RVT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( TRVT                );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( TRVT::ResultType    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( TRVT::TransposeType );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT::ResultType      );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT::TransposeType   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TVT::ResultType     );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TVT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RVT::ResultType     );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RVT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TRVT::ResultType    );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TRVT::TransposeType );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT::ElementType,   VT::ResultType::ElementType      );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( VT::ElementType,   VT::TransposeType::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( TVT::ElementType,  TVT::ResultType::ElementType     );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( TVT::ElementType,  TVT::TransposeType::ElementType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RVT::ElementType,  RVT::ResultType::ElementType     );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RVT::ElementType,  RVT::TransposeType::ElementType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( TRVT::ElementType, TRVT::ResultType::ElementType    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( TRVT::ElementType, TRVT::TransposeType::ElementType );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the alignment of different HybridVector instances.
//
// \return void
// \param type The string representation of the given template type.
// \exception std::runtime_error Error detected.
//
// This function performs a test of the alignment of a HybridVector instance of the given
// element type. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void ClassTest::testAlignment( const std::string& type )
{
   constexpr size_t SIMDSIZE ( blaze::SIMDTrait<Type>::size );
   constexpr size_t alignment( blaze::max( blaze::AlignmentOf_v<Type>, alignof( size_t ) ) );
   constexpr size_t overhead ( blaze::max( alignof(Type), alignof(size_t) ) );


   //=====================================================================================
   // Single vector alignment test (aligned/padded)
   //=====================================================================================

   {
      using AlignedPadded =
         blaze::HybridVector<Type,7UL,blaze::rowVector,blaze::aligned,blaze::padded>;

      BLAZE_STATIC_ASSERT( blaze::IsAligned_v<AlignedPadded> );
      BLAZE_STATIC_ASSERT( blaze::IsPadded_v<AlignedPadded> );
      BLAZE_STATIC_ASSERT( sizeof(AlignedPadded) == sizeof(Type)*blaze::nextMultiple( 7UL, SIMDSIZE ) + alignment );

      const AlignedPadded vec( 7UL );

      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );

      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: Vector alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Element type      : " << type << "\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation         : " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Single vector alignment test (aligned/unpadded)
   //=====================================================================================

   {
      using AlignedUnpadded =
         blaze::HybridVector<Type,7UL,blaze::rowVector,blaze::aligned,blaze::unpadded>;

      BLAZE_STATIC_ASSERT( blaze::IsAligned_v<AlignedUnpadded> );
      BLAZE_STATIC_ASSERT( !blaze::IsPadded_v<AlignedUnpadded> );
      BLAZE_STATIC_ASSERT( sizeof(AlignedUnpadded) == sizeof(Type)*blaze::nextMultiple( 7UL, SIMDSIZE ) + alignment );

      const AlignedUnpadded vec( 7UL );

      const size_t deviation( reinterpret_cast<size_t>( &vec[0] ) % alignment );

      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: Vector alignment test\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Element type      : " << type << "\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation         : " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Single vector alignment test (unaligned/padded)
   //=====================================================================================

   {
      using UnalignedPadded =
         blaze::HybridVector<Type,7UL,blaze::rowVector,blaze::unaligned,blaze::padded>;

      BLAZE_STATIC_ASSERT( !blaze::IsAligned_v<UnalignedPadded> );
      BLAZE_STATIC_ASSERT( blaze::IsPadded_v<UnalignedPadded> );
      BLAZE_STATIC_ASSERT( sizeof(UnalignedPadded) == sizeof(Type)*blaze::nextMultiple( 7UL, SIMDSIZE ) + overhead );
   }


   //=====================================================================================
   // Single vector alignment test (unaligned/unpadded)
   //=====================================================================================

   {
      using UnalignedUnpadded =
         blaze::HybridVector<Type,7UL,blaze::rowVector,blaze::unaligned,blaze::unpadded>;

      BLAZE_STATIC_ASSERT( !blaze::IsAligned_v<UnalignedUnpadded> );
      BLAZE_STATIC_ASSERT( !blaze::IsPadded_v<UnalignedUnpadded> );
      BLAZE_STATIC_ASSERT( sizeof(UnalignedUnpadded) == blaze::nextMultiple( sizeof(Type)*7UL, overhead ) + overhead );
   }


   //=====================================================================================
   // Static array alignment test (aligned/padded)
   //=====================================================================================

   {
      using AlignedPadded =
         blaze::HybridVector<Type,7UL,blaze::rowVector,blaze::aligned,blaze::padded>;

      const AlignedPadded init( 7UL );
      const std::array<AlignedPadded,7UL> vecs{ init, init, init, init, init, init, init };

      for( size_t i=0; i<vecs.size(); ++i )
      {
         const size_t deviation( reinterpret_cast<size_t>( &vecs[i][0] ) % alignment );

         if( deviation != 0UL ) {
            std::ostringstream oss;
            oss << " Test: Static array alignment test\n"
                << " Error: Invalid alignment at index " << i << " detected\n"
                << " Details:\n"
                << "   Element type      : " << type << "\n"
                << "   Expected alignment: " << alignment << "\n"
                << "   Deviation         : " << deviation << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Static array alignment test (aligned/unpadded)
   //=====================================================================================

   {
      using AlignedUnpadded =
         blaze::HybridVector<Type,7UL,blaze::rowVector,blaze::aligned,blaze::unpadded>;

      const AlignedUnpadded init( 7UL );
      const std::array<AlignedUnpadded,7UL> vecs{ init, init, init, init, init, init, init };

      for( size_t i=0; i<vecs.size(); ++i )
      {
         const size_t deviation( reinterpret_cast<size_t>( &vecs[i][0] ) % alignment );

         if( deviation != 0UL ) {
            std::ostringstream oss;
            oss << " Test: Static array alignment test\n"
                << " Error: Invalid alignment at index " << i << " detected\n"
                << " Details:\n"
                << "   Element type      : " << type << "\n"
                << "   Expected alignment: " << alignment << "\n"
                << "   Deviation         : " << deviation << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Dynamic array alignment test (aligned/padded)
   //=====================================================================================

   {
      using AlignedPadded =
         blaze::HybridVector<Type,7UL,blaze::rowVector,blaze::aligned,blaze::padded>;
      using AllocatorType = blaze::AlignedAllocator<AlignedPadded>;

      const AlignedPadded init( 7UL );
      const std::vector<AlignedPadded,AllocatorType> vecs( 7UL, init );

      for( size_t i=0; i<vecs.size(); ++i )
      {
         const size_t deviation( reinterpret_cast<size_t>( &vecs[i][0] ) % alignment );

         if( deviation != 0UL ) {
            std::ostringstream oss;
            oss << " Test: Dynamic array alignment test\n"
                << " Error: Invalid alignment at index " << i << " detected\n"
                << " Details:\n"
                << "   Element type      : " << type << "\n"
                << "   Expected alignment: " << alignment << "\n"
                << "   Deviation         : " << deviation << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Dynamic array alignment test (aligned/unpadded)
   //=====================================================================================

   {
      using AlignedUnpadded =
         blaze::HybridVector<Type,7UL,blaze::rowVector,blaze::aligned,blaze::unpadded>;
      using AllocatorType = blaze::AlignedAllocator<AlignedUnpadded>;

      const AlignedUnpadded init( 7UL );
      const std::vector<AlignedUnpadded,AllocatorType> vecs( 7UL, init );

      for( size_t i=0; i<vecs.size(); ++i )
      {
         const size_t deviation( reinterpret_cast<size_t>( &vecs[i][0] ) % alignment );

         if( deviation != 0UL ) {
            std::ostringstream oss;
            oss << " Test: Dynamic array alignment test\n"
                << " Error: Invalid alignment at index " << i << " detected\n"
                << " Details:\n"
                << "   Element type      : " << type << "\n"
                << "   Expected alignment: " << alignment << "\n"
                << "   Deviation         : " << deviation << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the size of the given hybrid vector.
//
// \param vector The hybrid vector to be checked.
// \param expectedSize The expected size of the hybrid vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the size of the given hybrid vector. In case the actual size
// does not correspond to the given expected size, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >  // Type of the hybrid vector
void ClassTest::checkSize( const Type& vector, size_t expectedSize ) const
{
   if( size( vector ) != expectedSize ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid size detected\n"
          << " Details:\n"
          << "   Size         : " << size( vector ) << "\n"
          << "   Expected size: " << expectedSize << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the capacity of the given hybrid vector.
//
// \param vector The hybrid vector to be checked.
// \param minCapacity The expected minimum capacity of the hybrid vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of the given hybrid vector. In case the actual capacity
// is smaller than the given expected minimum capacity, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >  // Type of the hybrid vector
void ClassTest::checkCapacity( const Type& vector, size_t minCapacity ) const
{
   if( capacity( vector ) < minCapacity ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Capacity                 : " << capacity( vector ) << "\n"
          << "   Expected minimum capacity: " << minCapacity << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of non-zero elements of the given hybrid vector.
//
// \param vector The hybrid vector to be checked.
// \param expectedNonZeros The expected number of non-zero elements of the hybrid vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements of the given hybrid vector. In
// case the actual number of non-zero elements does not correspond to the given expected
// number, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the hybrid vector
void ClassTest::checkNonZeros( const Type& vector, size_t expectedNonZeros ) const
{
   if( nonZeros( vector ) != expectedNonZeros ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of non-zero elements\n"
          << " Details:\n"
          << "   Number of non-zeros         : " << nonZeros( vector ) << "\n"
          << "   Expected number of non-zeros: " << expectedNonZeros << "\n";
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
/*!\brief Testing the functionality of the HybridVector class template.
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
/*!\brief Macro for the execution of the HybridVector class test.
*/
#define RUN_HYBRIDVECTOR_CLASS_TEST \
   blazetest::mathtest::vectors::hybridvector::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace hybridvector

} // namespace vectors

} // namespace mathtest

} // namespace blazetest

#endif
