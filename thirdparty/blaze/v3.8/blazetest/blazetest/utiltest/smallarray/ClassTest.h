//=================================================================================================
/*!
//  \file blazetest/utiltest/smallarray/ClassTest.h
//  \brief Header file for the SmallArray test
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

#ifndef _BLAZETEST_UTILTEST_SMALLARRAY_CLASSTEST_H_
#define _BLAZETEST_UTILTEST_SMALLARRAY_CLASSTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <list>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/Random.h>
#include <blaze/util/SmallArray.h>
#include <blazetest/utiltest/IntResource.h>


namespace blazetest {

namespace utiltest {

namespace smallarray {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for the test of the SmallArray class template.
//
// This class represents the collection of tests for the SmallArray class template.
*/
template< size_t N >  // Number of preallocated elements
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
   void testConstructors();
   void testAssignment  ();
   void testSubscript   ();
   void testAt          ();
   void testIterator    ();
   void testClear       ();
   void testResize      ();
   void testReserve     ();
   void testShrinkToFit ();
   void testPushBack    ();
   void testInsert      ();
   void testErase       ();
   void testSwap        ();

   template< typename Type >
   void checkSize( const Type& array, size_t expectedSize ) const;

   template< typename Type >
   void checkCapacity( const Type& array, size_t minCapacity ) const;

   void checkCount( size_t expectedCount ) const;
   //@}
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   using InitList = std::initializer_list<IntResource>;  //!< Initializer list of integer resources.
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
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global output operator for small arrays.
//
// \param os Reference to the output stream.
// \param v Reference to a constant small array object.
// \return Reference to the output stream.
*/
template< typename Type  // Data type of the elements
        , size_t N >     // Number of preallocated elements
inline std::ostream& operator<<( std::ostream& os, const blaze::SmallArray<Type,N>& sv )
{
   os << "(";
   for( const Type& value : sv )
      os << " " << value;
   os << " )";

   return os;
}
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the SmallArray class test.
//
// \exception std::runtime_error Operation error detected.
*/
template< size_t N >  // Number of preallocated elements
ClassTest<N>::ClassTest()
{
   testConstructors();
   testAssignment();
   testSubscript();
   testAt();
   testIterator();
   testClear();
   testResize();
   testReserve();
   testShrinkToFit();
   testPushBack();
   testInsert();
   testErase();
   testSwap();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the SmallArray constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the SmallArray class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testConstructors()
{
   //=====================================================================================
   // Default constructor
   //=====================================================================================

   {
      test_ = "SmallArray default constructor";

      blaze::SmallArray<IntResource,N> arr;

      checkSize ( arr, 0UL );
      checkCount( 0U );
   }


   //=====================================================================================
   // Size constructor
   //=====================================================================================

   {
      test_ = "SmallArray size constructor (size 0)";

      blaze::SmallArray<IntResource,N> arr( 0UL );

      checkSize ( arr, 0UL );
      checkCount( 0U );
   }

   {
      test_ = "SmallArray size constructor (size 4)";

      blaze::SmallArray<IntResource,N> arr( 4UL );

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );
      checkCount   ( 4U );
   }

   {
      test_ = "SmallArray size constructor (size 5)";

      blaze::SmallArray<IntResource,N> arr( 5UL );

      checkSize    ( arr, 5UL );
      checkCapacity( arr, 5UL );
      checkCount   ( 5U );
   }

   {
      test_ = "SmallArray size constructor (size 6)";

      blaze::SmallArray<IntResource,N> arr( 6UL );

      checkSize    ( arr, 6UL );
      checkCapacity( arr, 6UL );
      checkCount   ( 6U );
   }


   //=====================================================================================
   // Homogeneous initialization
   //=====================================================================================

   {
      test_ = "SmallArray homogeneous initialization constructor (size 0)";

      blaze::SmallArray<IntResource,N> arr( 0UL, 2 );

      checkSize ( arr, 0UL );
      checkCount( 0U );
   }

   {
      test_ = "SmallArray homogeneous initialization constructor (size 4)";

      blaze::SmallArray<IntResource,N> arr( 4UL, 2 );

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );
      checkCount   ( 4U );

      if( arr[0] != 2 || arr[1] != 2 || arr[2] != 2 || arr[3] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray homogeneous initialization constructor (size 5)";

      blaze::SmallArray<IntResource,N> arr( 5UL, 2 );

      checkSize    ( arr, 5UL );
      checkCapacity( arr, 5UL );
      checkCount   ( 5U );

      if( arr[0] != 2 || arr[1] != 2 || arr[2] != 2 || arr[3] != 2 || arr[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray homogeneous initialization constructor (size 6)";

      blaze::SmallArray<IntResource,N> arr( 6UL, 2 );

      checkSize    ( arr, 6UL );
      checkCapacity( arr, 6UL );
      checkCount   ( 6U );

      if( arr[0] != 2 || arr[1] != 2 || arr[2] != 2 || arr[3] != 2 || arr[4] != 2 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 2 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Range initialization
   //=====================================================================================

   {
      test_ = "SmallArray range constructor (size 4)";

      std::list<int> list{ 1, 2, 3, 4 };
      blaze::SmallArray<IntResource,N> arr( list.begin(), list.end() );

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );
      checkCount   ( 4U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray range constructor (size 5)";

      std::list<int> list{ 1, 2, 3, 4, 5 };
      blaze::SmallArray<IntResource,N> arr( list.begin(), list.end() );

      checkSize    ( arr, 5UL );
      checkCapacity( arr, 5UL );
      checkCount   ( 5U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray range constructor (size 6)";

      std::list<int> list{ 1, 2, 3, 4, 5, 6 };
      blaze::SmallArray<IntResource,N> arr( list.begin(), list.end() );

      checkSize    ( arr, 6UL );
      checkCapacity( arr, 6UL );
      checkCount   ( 6U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 || arr[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // List initialization
   //=====================================================================================

   {
      test_ = "SmallArray initializer list constructor (size 4)";

      blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4 } );

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );
      checkCount   ( 4U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray initializer list constructor (size 5)";

      blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5 } );

      checkSize    ( arr, 5UL );
      checkCapacity( arr, 5UL );
      checkCount   ( 5U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray initializer list constructor (size 6)";

      blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5, 6 } );

      checkSize    ( arr, 6UL );
      checkCapacity( arr, 6UL );
      checkCount   ( 6U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 || arr[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy constructor
   //=====================================================================================

   {
      test_ = "SmallArray copy constructor (size 0)";

      blaze::SmallArray<IntResource,N> arr1( 0UL );
      blaze::SmallArray<IntResource,N> arr2( arr1 );

      checkSize    ( arr2, 0UL );
      checkCapacity( arr2, 0UL );
      checkCount   ( 0U );
   }

   {
      test_ = "SmallArray copy constructor (size 4)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4 } );
      blaze::SmallArray<IntResource,N> arr2( arr1 );

      checkSize    ( arr2, 4UL );
      checkCapacity( arr2, 4UL );
      checkCount   ( 8U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray copy constructor (size 5)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4, 5 } );
      blaze::SmallArray<IntResource,N> arr2( arr1 );

      checkSize    ( arr2, 5UL );
      checkCapacity( arr2, 5UL );
      checkCount   ( 10U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 || arr2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray copy constructor (size 6)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4, 5, 6 } );
      blaze::SmallArray<IntResource,N> arr2( arr1 );

      checkSize    ( arr2, 6UL );
      checkCapacity( arr2, 6UL );
      checkCount   ( 12U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 || arr2[4] != 5 || arr2[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Move constructor
   //=====================================================================================

   {
      test_ = "SmallArray move constructor (size 0)";

      blaze::SmallArray<IntResource,N> arr1( 0UL );
      blaze::SmallArray<IntResource,N> arr2( std::move( arr1 ) );

      checkSize    ( arr2, 0UL );
      checkCapacity( arr2, 0UL );
      checkCount   ( 0U );
   }

   {
      test_ = "SmallArray move constructor (size 4)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4 } );
      blaze::SmallArray<IntResource,N> arr2( std::move( arr1 ) );

      checkSize    ( arr2, 4UL );
      checkCapacity( arr2, 4UL );
      checkCount   ( 4U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray move constructor (size 5)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4, 5 } );
      blaze::SmallArray<IntResource,N> arr2( std::move( arr1 ) );

      checkSize    ( arr2, 5UL );
      checkCapacity( arr2, 5UL );
      checkCount   ( 5U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 || arr2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray move constructor (size 6)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4, 5, 6 } );
      blaze::SmallArray<IntResource,N> arr2( std::move( arr1 ) );

      checkSize    ( arr2, 6UL );
      checkCapacity( arr2, 6UL );
      checkCount   ( 6U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 || arr2[4] != 5 || arr2[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SmallArray assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the SmallArray class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testAssignment()
{
   //=====================================================================================
   // List assignment
   //=====================================================================================

   {
      test_ = "SmallArray initializer list assignment (size 3 to 4)";

      blaze::SmallArray<IntResource,N> arr( InitList{ 11, 12, 13 } );
      arr = InitList{ 1, 2, 3, 4 };

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );
      checkCount   ( 4U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray initializer list assignment (size 8 to 4)";

      blaze::SmallArray<IntResource,N> arr( InitList{ 11, 12, 13, 14, 15, 16, 17, 18 } );
      arr = InitList{ 1, 2, 3, 4 };

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );
      checkCount   ( 4U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray initializer list assignment (size 3 to 5)";

      blaze::SmallArray<IntResource,N> arr( InitList{ 11, 12, 13 } );
      arr = InitList{ 1, 2, 3, 4, 5 };

      checkSize    ( arr, 5UL );
      checkCapacity( arr, 5UL );
      checkCount   ( 5U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray initializer list assignment (size 8 to 5)";

      blaze::SmallArray<IntResource,N> arr( InitList{ 11, 12, 13, 14, 15, 16, 17, 18 } );
      arr = InitList{ 1, 2, 3, 4, 5 };

      checkSize    ( arr, 5UL );
      checkCapacity( arr, 5UL );
      checkCount   ( 5U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray initializer list assignment (size 3 to 6)";

      blaze::SmallArray<IntResource,N> arr( InitList{ 11, 12, 13 } );
      arr = InitList{ 1, 2, 3, 4, 5, 6 };

      checkSize    ( arr, 6UL );
      checkCapacity( arr, 6UL );
      checkCount   ( 6U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 || arr[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray initializer list assignment (size 8 to 6)";

      blaze::SmallArray<IntResource,N> arr( InitList{ 11, 12, 13, 14, 15, 16, 17, 18 } );
      arr = InitList{ 1, 2, 3, 4, 5, 6 };

      checkSize    ( arr, 6UL );
      checkCapacity( arr, 6UL );
      checkCount   ( 6U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 || arr[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "SmallArray copy assignment (size 4)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4 } );
      blaze::SmallArray<IntResource,N> arr2;
      arr2 = arr1;

      checkSize    ( arr2, 4UL );
      checkCapacity( arr2, 4UL );
      checkCount   ( 8U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray copy assignment (size 5)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4, 5 } );
      blaze::SmallArray<IntResource,N> arr2;
      arr2 = arr1;

      checkSize    ( arr2, 5UL );
      checkCapacity( arr2, 5UL );
      checkCount   ( 10U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 || arr2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray copy assignment (size 6)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4, 5, 6 } );
      blaze::SmallArray<IntResource,N> arr2;
      arr2 = arr1;

      checkSize    ( arr2, 6UL );
      checkCapacity( arr2, 6UL );
      checkCount   ( 12U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 || arr2[4] != 5 || arr2[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray copy assignment stress test";

      blaze::SmallArray<int,N> arr1;
      const int min( -10 );
      const int max(  10 );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t size( blaze::rand<size_t>( 0UL, 10UL ) );
         blaze::SmallArray<int,N> arr2( size );
         for( int& element : arr2 ) {
            element = blaze::rand<int>( min, max );
         }

         arr1 = arr2;

         if( arr1 != arr2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << arr1 << "\n"
                << "   Expected result:\n" << arr2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Move assignment
   //=====================================================================================

   {
      test_ = "SmallArray move assignment (size 3 to 4)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4 } );
      blaze::SmallArray<IntResource,N> arr2( InitList{ 11, 12, 13 } );

      arr2 = std::move( arr1 );

      checkSize    ( arr2, 4UL );
      checkCapacity( arr2, 4UL );
      checkCount   ( 4U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray move assignment (size 8 to 4)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4 } );
      blaze::SmallArray<IntResource,N> arr2( InitList{ 11, 12, 13, 14, 15, 16, 17, 18 } );

      arr2 = std::move( arr1 );

      checkSize    ( arr2, 4UL );
      checkCapacity( arr2, 4UL );
      checkCount   ( 4U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray move assignment (size 3 to 5)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4, 5 } );
      blaze::SmallArray<IntResource,N> arr2( InitList{ 11, 12, 13 } );

      arr2 = std::move( arr1 );

      checkSize    ( arr2, 5UL );
      checkCapacity( arr2, 5UL );
      checkCount   ( 5U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 || arr2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray move assignment (size 8 to 5)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4, 5 } );
      blaze::SmallArray<IntResource,N> arr2( InitList{ 11, 12, 13, 14, 15, 16, 17, 18 } );

      arr2 = std::move( arr1 );

      checkSize    ( arr2, 5UL );
      checkCapacity( arr2, 5UL );
      checkCount   ( 5U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 || arr2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray move assignment (size 3 to 6)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4, 5, 6 } );
      blaze::SmallArray<IntResource,N> arr2( InitList{ 11, 12, 13 } );

      arr2 = std::move( arr1 );

      checkSize    ( arr2, 6UL );
      checkCapacity( arr2, 6UL );
      checkCount   ( 6U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 || arr2[4] != 5 || arr2[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray move assignment (size 8 to 6)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4, 5, 6 } );
      blaze::SmallArray<IntResource,N> arr2( InitList{ 11, 12, 13, 14, 15, 16, 17, 18 } );

      arr2 = std::move( arr1 );

      checkSize    ( arr2, 6UL );
      checkCapacity( arr2, 6UL );
      checkCount   ( 6U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 || arr2[4] != 5 || arr2[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << arr2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SmallArray subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the SmallArray class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testSubscript()
{
   {
      test_ = "SmallArray::operator[] (size 4)";

      // Assignment to the element at index 2
      blaze::SmallArray<int,N> arr{ 0, 0, 1, 0 };

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[2] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 0 0 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 3
      arr[3] = 3;

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[2] != 1 || arr[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 0 0 1 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 0
      arr[0] = 4;

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[0] != 4 || arr[2] != 1 || arr[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 0 1 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element at index 2
      arr[2] += arr[3];

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[0] != 4 || arr[2] != 4 || arr[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 0 4 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element at index 1
      arr[1] -= 2;

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[0] != 4 || arr[1] != -2 || arr[2] != 4 || arr[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 -2 4 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element at index 3
      arr[3] *= -3;

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[0] != 4 || arr[1] != -2 || arr[2] != 4 || arr[3] != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 -2 4 -9 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element at index 2
      arr[2] /= 2;

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[0] != 4 || arr[1] != -2 || arr[2] != 2 || arr[3] != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 -2 2 -9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray::operator[] (size 7)";

      // Assignment to the element at index 2
      blaze::SmallArray<int,N> arr{ 0, 0, 1, 0, 0, 0, 0 };

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[2] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 0 0 1 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 5
      arr[5] = 2;

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[2] != 1 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 0 0 1 0 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 3
      arr[3] = 3;

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[2] != 1 || arr[3] != 3 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 0 0 1 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 0
      arr[0] = 4;

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[0] != 4 || arr[2] != 1 || arr[3] != 3 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 0 1 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element at index 2
      arr[2] += arr[3];

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[0] != 4 || arr[2] != 4 || arr[3] != 3 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 0 4 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element at index 1
      arr[1] -= arr[5];

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[0] != 4 || arr[1] != -2 || arr[2] != 4 || arr[3] != 3 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 -2 4 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element at index 3
      arr[3] *= -3;

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[0] != 4 || arr[1] != -2 || arr[2] != 4 || arr[3] != -9 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 -2 4 -9 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element at index 2
      arr[2] /= 2;

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[0] != 4 || arr[1] != -2 || arr[2] != 2 || arr[3] != -9 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 -2 2 -9 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the SmallArray class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the SmallArray class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testAt()
{
   {
      test_ = "SmallArray::at() (size 4)";

      // Assignment to the element at index 2
      blaze::SmallArray<int,N> arr{ 0, 0, 1, 0 };

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr.at(2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 0 0 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 3
      arr.at(3) = 3;

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[2] != 1 || arr[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 0 0 1 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 0
      arr.at(0) = 4;

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[0] != 4 || arr[2] != 1 || arr[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 0 1 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element at index 2
      arr.at(2) += arr.at(3);

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[0] != 4 || arr[2] != 4 || arr[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 0 4 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element at index 1
      arr.at(1) -= 2;

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[0] != 4 || arr[1] != -2 || arr[2] != 4 || arr[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 -2 4 3 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element at index 3
      arr.at(3) *= -3;

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[0] != 4 || arr[1] != -2 || arr[2] != 4 || arr[3] != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 -2 4 -9 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element at index 2
      arr.at(2) /= 2;

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );

      if( arr[0] != 4 || arr[1] != -2 || arr[2] != 2 || arr[3] != -9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 -2 2 -9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray::at() (size 7)";

      // Assignment to the element at index 2
      blaze::SmallArray<int,N> arr{ 0, 0, 1, 0, 0, 0, 0 };

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[2] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 0 0 1 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 5
      arr.at(5) = 2;

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[2] != 1 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 0 0 1 0 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 3
      arr.at(3) = 3;

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[2] != 1 || arr[3] != 3 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 0 0 1 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Assignment to the element at index 0
      arr.at(0) = 4;

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[0] != 4 || arr[2] != 1 || arr[3] != 3 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 0 1 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Addition assignment to the element at index 2
      arr.at(2) += arr.at(3);

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[0] != 4 || arr[2] != 4 || arr[3] != 3 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 0 4 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Subtraction assignment to the element at index 1
      arr.at(1) -= arr.at(5);

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[0] != 4 || arr[1] != -2 || arr[2] != 4 || arr[3] != 3 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 -2 4 3 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Multiplication assignment to the element at index 3
      arr.at(3) *= -3;

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[0] != 4 || arr[1] != -2 || arr[2] != 4 || arr[3] != -9 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 -2 4 -9 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Division assignment to the element at index 2
      arr.at(2) /= 2;

      checkSize    ( arr, 7UL );
      checkCapacity( arr, 7UL );

      if( arr[0] != 4 || arr[1] != -2 || arr[2] != 2 || arr[3] != -9 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 4 -2 2 -9 0 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SmallArray iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the SmallArray class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testIterator()
{
   using ArrayType     = blaze::SmallArray<int,N>;
   using Iterator      = typename ArrayType::Iterator;
   using ConstIterator = typename ArrayType::ConstIterator;

   ArrayType arr{ 1, 0, -2, -3 };

   // Testing the Iterator default constructor
   {
      test_ = "Iterator default constructor";

      Iterator it{};

      if( it != Iterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing the ConstIterator default constructor
   {
      test_ = "ConstIterator default constructor";

      ConstIterator it{};

      if( it != ConstIterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing conversion from Iterator to ConstIterator
   {
      test_ = "Iterator/ConstIterator conversion";

      ConstIterator it( begin( arr ) );

      if( it == end( arr ) || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator conversion detected\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements via Iterator
   {
      test_ = "Iterator subtraction";

      const size_t number( end( arr ) - begin( arr ) );

      if( number != 4UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements via ConstIterator
   {
      test_ = "ConstIterator subtraction";

      const size_t number( cend( arr ) - cbegin( arr ) );

      if( number != 4UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing read-only access via ConstIterator
   {
      test_ = "Read-only access via ConstIterator";

      ConstIterator it ( cbegin( arr ) );
      ConstIterator end( cend( arr ) );

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid initial iterator detected\n";
         throw std::runtime_error( oss.str() );
      }

      ++it;

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      --it;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it++;

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it--;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it += 2UL;

      if( it == end || *it != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator addition assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it -= 2UL;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator subtraction assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it + 3UL;

      if( it == end || *it != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar addition failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it - 3UL;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar subtraction failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = 4UL + it;

      if( it != end ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scalar/iterator addition failed\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing assignment via Iterator
   {
      test_ = "Assignment via Iterator";

      int value = 6;

      for( Iterator it=begin( arr ); it!=end( arr ); ++it ) {
         *it = value++;
      }

      if( arr[0] != 6 || arr[1] != 7 || arr[2] != 8 || arr[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 6 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing addition assignment via Iterator
   {
      test_ = "Addition assignment via Iterator";

      int value = 2;

      for( Iterator it=begin( arr ); it!=end( arr ); ++it ) {
         *it += value++;
      }

      if( arr[0] != 8 || arr[1] != 10 || arr[2] != 12 || arr[3] != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 8 10 12 14 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing subtraction assignment via Iterator
   {
      test_ = "Subtraction assignment via Iterator";

      int value = 2;

      for( Iterator it=begin( arr ); it!=end( arr ); ++it ) {
         *it -= value++;
      }

      if( arr[0] != 6 || arr[1] != 7 || arr[2] != 8 || arr[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 6 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing multiplication assignment via Iterator
   {
      test_ = "Multiplication assignment via Iterator";

      int value = 1;

      for( Iterator it=begin( arr ); it!=end( arr ); ++it ) {
         *it *= value++;
      }

      if( arr[0] != 6 || arr[1] != 14 || arr[2] != 24 || arr[3] != 36 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 6 14 24 36 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing division assignment via Iterator
   {
      test_ = "Division assignment via Iterator";

      for( Iterator it=begin( arr ); it!=end( arr ); ++it ) {
         *it /= 2;
      }

      if( arr[0] != 3 || arr[1] != 7 || arr[2] != 12 || arr[3] != 18 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 3 7 12 18 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() member function of the SmallArray class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() member function of the SmallArray class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testClear()
{
   test_ = "SmallArray::clear()";

   // Clearing a default constructed array
   {
      blaze::SmallArray<IntResource,N> arr;

      clear( arr );

      checkSize ( arr, 0UL );
      checkCount( 0U );
   }

   // Clearing an initialized array
   {
      // Initialization check
      blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4 } );

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );
      checkCount   ( 4U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the array
      clear( arr );

      checkSize ( arr, 0UL );
      checkCount( 0U );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c resize() member function of the SmallArray class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c resize() member function of the SmallArray class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testResize()
{
   {
      test_ = "SmallArray::resize( size_t )";

      // Initialization check
      blaze::SmallArray<IntResource,N> arr;

      checkSize ( arr, 0UL );
      checkCount( 0U );

      // Resizing to 0
      arr.resize( 0UL );

      checkSize ( arr, 0UL );
      checkCount( 0U );

      // Resizing to 4
      arr.resize( 4UL );
      arr[0] = 1;
      arr[1] = 2;
      arr[2] = 3;
      arr[3] = 4;

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );
      checkCount   ( 4U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 6
      arr.resize( 6UL );
      arr[4] = 5;
      arr[5] = 6;

      checkSize    ( arr, 6UL );
      checkCapacity( arr, 6UL );
      checkCount   ( 6U );

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 || arr[5] != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 3
      arr.resize( 3UL );
      arr[0] = 11;
      arr[1] = 12;
      arr[2] = 13;

      checkSize    ( arr, 3UL );
      checkCapacity( arr, 3UL );
      checkCount   ( 3U );

      if( arr[0] != 11 || arr[1] != 12 || arr[2] != 13 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 11 12 13 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 6
      arr.resize( 6UL );
      arr[3] = 14;
      arr[4] = 15;
      arr[5] = 16;

      checkSize    ( arr, 6UL );
      checkCapacity( arr, 6UL );
      checkCount   ( 6U );

      if( arr[0] != 11 || arr[1] != 12 || arr[2] != 13 || arr[3] != 14 || arr[4] != 15 || arr[5] != 16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 11 12 13 14 15 16 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0
      arr.resize( 0UL );

      checkSize ( arr, 0UL );
      checkCount( 0U );
   }

   {
      test_ = "SmallArray::resize( size_t, const Type& )";

      // Initialization check
      blaze::SmallArray<IntResource,N> arr;

      checkSize ( arr, 0UL );
      checkCount( 0U );

      // Resizing to 0
      arr.resize( 0UL, 2 );

      checkSize ( arr, 0UL );
      checkCount( 0U );

      // Resizing to 4
      arr.resize( 4UL, 2 );

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );
      checkCount   ( 4U );

      if( arr[0] != 2 || arr[1] != 2 || arr[2] != 2 || arr[3] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 6
      arr.resize( 6UL, 2 );

      checkSize    ( arr, 6UL );
      checkCapacity( arr, 6UL );
      checkCount   ( 6U );

      if( arr[0] != 2 || arr[1] != 2 || arr[2] != 2 || arr[3] != 2 || arr[4] != 2 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 2 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 3
      arr.resize( 3UL, 2 );

      checkSize    ( arr, 3UL );
      checkCapacity( arr, 3UL );
      checkCount   ( 3U );

      if( arr[0] != 2 || arr[1] != 2 || arr[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 6
      arr.resize( 6UL, 2 );

      checkSize    ( arr, 6UL );
      checkCapacity( arr, 6UL );
      checkCount   ( 6U );

      if( arr[0] != 2 || arr[1] != 2 || arr[2] != 2 || arr[3] != 2 || arr[4] != 2 || arr[5] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Resize operation failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 2 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resizing to 0
      arr.resize( 0UL, 2 );

      checkSize ( arr, 0UL );
      checkCount( 0U );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the SmallArray class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the SmallArray class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testReserve()
{
   test_ = "SmallArray::reserve()";

   // Initialization check
   blaze::SmallArray<IntResource,N> arr;

   checkSize ( arr, 0UL );
   checkCount( 0U );

   // Increasing the capacity of the array
   arr.reserve( 4UL );

   checkSize    ( arr, 0UL );
   checkCapacity( arr, 4UL );
   checkCount   ( 0U );

   // Further increasing the capacity of the array
   arr.reserve( 8UL );

   checkSize    ( arr, 0UL );
   checkCapacity( arr, 8UL );
   checkCount   ( 0U );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c shrinkToFit() member function of the SmallArray class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c shrinkToFit() member function of the SmallArray class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testShrinkToFit()
{
   test_ = "SmallArray::shrinkToFit()";

   // Shrinking an array without excessive capacity
   {
      blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4 } );

      arr.shrinkToFit();

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );
      checkCount   ( 4U );

      if( arr.capacity() > blaze::max( 4UL, N ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the array failed\n"
             << " Details:\n"
             << "   Capacity: " << arr.capacity() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the array failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Shrinking an array with excessive capacity (size 4)
   {
      blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4 } );
      arr.reserve( 100UL );

      arr.shrinkToFit();

      checkSize    ( arr, 4UL );
      checkCapacity( arr, 4UL );
      checkCount   ( 4U );

      if( arr.capacity() > blaze::max( 4UL, N ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the array failed\n"
             << " Details:\n"
             << "   Capacity: " << arr.capacity() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the array failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Shrinking an array with excessive capacity (size 8)
   {
      blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5, 6, 7, 8 } );
      arr.reserve( 100UL );

      arr.shrinkToFit();

      checkSize    ( arr, 8UL );
      checkCapacity( arr, 8UL );
      checkCount   ( 8U );

      if( arr.capacity() > blaze::max( 8UL, N ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the array failed\n"
             << " Details:\n"
             << "   Capacity: " << arr.capacity() << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ||
          arr[4] != 5 || arr[5] != 6 || arr[6] != 7 || arr[7] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Shrinking the array failed\n"
             << " Details:\n"
             << "   Result:\n" << arr << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c pushBack() member function of the SmallArray class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c pushBack() member function of the SmallArray class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testPushBack()
{
   test_ = "SmallArray::pushBack() (size 4)";

   blaze::SmallArray<IntResource,N> arr;

   checkSize( arr, 0UL );

   arr.pushBack( 1 );
   arr.pushBack( 2 );
   arr.pushBack( 3 );
   arr.pushBack( 4 );
   arr.pushBack( 5 );

   checkSize    ( arr, 5UL );
   checkCapacity( arr, 5UL );
   checkCount   ( 5U );

   if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << arr << "\n"
          << "   Expected result:\n( 1 2 3 4 5 )\n";
      throw std::runtime_error( oss.str() );
   }

   arr.pushBack( 6 );
   arr.pushBack( 7 );
   arr.pushBack( 8 );

   checkSize    ( arr, 8UL );
   checkCapacity( arr, 8UL );
   checkCount   ( 8U );

   if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ||
       arr[4] != 5 || arr[5] != 6 || arr[6] != 7 || arr[7] != 8 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << arr << "\n"
          << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the SmallArray class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the SmallArray class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testInsert()
{
   {
      // Inserting into an empty small array
      {
         test_ = "SmallArray::insert( Iterator, const Type& ) (empty array)";

         blaze::SmallArray<IntResource,N> arr;
         int value = 1;

         auto pos = arr.insert( arr.begin(), value );

         checkSize    ( arr, 1UL );
         checkCapacity( arr, 1UL );
         checkCount   ( 1U );

         if( pos == arr.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the beginning of a small array (x 2 3 4)
      {
         test_ = "SmallArray::insert( Iterator, const Type& ) (x 2 3 4)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 2, 3, 4 } );
         int value = 1;

         auto pos = arr.insert( arr.begin(), value );

         checkSize    ( arr, 4UL );
         checkCapacity( arr, 4UL );
         checkCount   ( 4U );

         if( pos == arr.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting in the middle of a small array (1 x 3 4)
      {
         test_ = "SmallArray::insert( Iterator, const Type& ) (1 x 3 4)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 3, 4 } );
         int value = 2;

         auto pos = arr.insert( arr.begin()+1UL, value );

         checkSize    ( arr, 4UL );
         checkCapacity( arr, 4UL );
         checkCount   ( 4U );

         if( pos == arr.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the end of a small array (1 2 3 x)
      {
         test_ = "SmallArray::insert( Iterator, const Type& ) (1 2 3 x)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3 } );
         int value = 4;

         auto pos = arr.insert( arr.end(), value );

         checkSize    ( arr, 4UL );
         checkCapacity( arr, 4UL );
         checkCount   ( 4U );

         if( pos == arr.end() || *pos != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the beginning of a small array (x 2 3 4 5 6)
      {
         test_ = "SmallArray::insert( Iterator, const Type& ) (x 2 3 4 5 6)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 2, 3, 4, 5, 6 } );
         int value = 1;

         auto pos = arr.insert( arr.begin(), value );

         checkSize    ( arr, 6UL );
         checkCapacity( arr, 6UL );
         checkCount   ( 6U );

         if( pos == arr.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 || arr[5] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting in the middle of a small array (1 x 3 4 5 6)
      {
         test_ = "SmallArray::insert( Iterator, const Type& ) (1 x 3 4 5 6)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 3, 4, 5, 6 } );
         int value = 2;

         auto pos = arr.insert( arr.begin()+1UL, value );

         checkSize    ( arr, 6UL );
         checkCapacity( arr, 6UL );
         checkCount   ( 6U );

         if( pos == arr.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 || arr[5] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the end of a small array (1 2 3 4 5 x)
      {
         test_ = "SmallArray::insert( Iterator, const Type& ) (1 2 3 4 5 x)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5 } );
         int value = 6;

         auto pos = arr.insert( arr.end(), value );

         checkSize    ( arr, 6UL );
         checkCapacity( arr, 6UL );
         checkCount   ( 6U );

         if( pos == arr.end() || *pos != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 6\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 || arr[5] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the beginning of a small array (x 2 3 4 5 6 7 8)
      {
         test_ = "SmallArray::insert( Iterator, const Type& ) (x 2 3 4 5 6 7 8)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 2, 3, 4, 5, 6, 7, 8 } );
         int value = 1;

         auto pos = arr.insert( arr.begin(), value );

         checkSize    ( arr, 8UL );
         checkCapacity( arr, 8UL );
         checkCount   ( 8U );

         if( pos == arr.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ||
             arr[4] != 5 || arr[5] != 6 || arr[6] != 7 || arr[7] != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting in the middle of a small array (1 x 3 4 5 6 7 8 )
      {
         test_ = "SmallArray::insert( Iterator, const Type& ) (1 x 3 4 5 6 7 8)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 3, 4, 5, 6, 7, 8 } );
         int value = 2;

         auto pos = arr.insert( arr.begin()+1UL, value );

         checkSize    ( arr, 8UL );
         checkCapacity( arr, 8UL );
         checkCount   ( 8U );

         if( pos == arr.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ||
             arr[4] != 5 || arr[5] != 6 || arr[6] != 7 || arr[7] != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the end of a small array (1 2 3 4 5 6 7 x)
      {
         test_ = "SmallArray::insert( Iterator, const Type& ) (1 2 3 4 5 6 7 x)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5, 6, 7 } );
         int value = 8;

         auto pos = arr.insert( arr.end(), value );

         checkSize    ( arr, 8UL );
         checkCapacity( arr, 8UL );
         checkCount   ( 8U );

         if( pos == arr.end() || *pos != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 8\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ||
             arr[4] != 5 || arr[5] != 6 || arr[6] != 7 || arr[7] != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      // Inserting into an empty small array
      {
         test_ = "SmallArray::insert( Iterator, Type&& ) (empty array)";

         blaze::SmallArray<IntResource,N> arr;

         auto pos = arr.insert( arr.begin(), 1 );

         checkSize    ( arr, 1UL );
         checkCapacity( arr, 1UL );
         checkCount   ( 1U );

         if( pos == arr.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the beginning of a small array (x 2 3 4)
      {
         test_ = "SmallArray::insert( Iterator, Type&& ) (x 2 3 4)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 2, 3, 4 } );

         auto pos = arr.insert( arr.begin(), 1 );

         checkSize    ( arr, 4UL );
         checkCapacity( arr, 4UL );
         checkCount   ( 4U );

         if( pos == arr.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting in the middle of a small array (1 x 3 4)
      {
         test_ = "SmallArray::insert( Iterator, Type&& ) (1 x 3 4)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 3, 4 } );

         auto pos = arr.insert( arr.begin()+1UL, 2 );

         checkSize    ( arr, 4UL );
         checkCapacity( arr, 4UL );
         checkCount   ( 4U );

         if( pos == arr.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the end of a small array (1 2 3 x)
      {
         test_ = "SmallArray::insert( Iterator, Type&& ) (1 2 3 x)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3 } );

         auto pos = arr.insert( arr.end(), 4 );

         checkSize    ( arr, 4UL );
         checkCapacity( arr, 4UL );
         checkCount   ( 4U );

         if( pos == arr.end() || *pos != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the beginning of a small array (x 2 3 4 5 6)
      {
         test_ = "SmallArray::insert( Iterator, Type&& ) (x 2 3 4 5 6)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 2, 3, 4, 5, 6 } );

         auto pos = arr.insert( arr.begin(), 1 );

         checkSize    ( arr, 6UL );
         checkCapacity( arr, 6UL );
         checkCount   ( 6U );

         if( pos == arr.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 || arr[5] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting in the middle of a small array (1 x 3 4 5 6)
      {
         test_ = "SmallArray::insert( Iterator, Type&& ) (1 x 3 4 5 6)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 3, 4, 5, 6 } );

         auto pos = arr.insert( arr.begin()+1UL, 2 );

         checkSize    ( arr, 6UL );
         checkCapacity( arr, 6UL );
         checkCount   ( 6U );

         if( pos == arr.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 || arr[5] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the end of a small array (1 2 3 4 5 x)
      {
         test_ = "SmallArray::insert( Iterator, Type&& ) (1 2 3 4 5 x)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5 } );

         auto pos = arr.insert( arr.end(), 6 );

         checkSize    ( arr, 6UL );
         checkCapacity( arr, 6UL );
         checkCount   ( 6U );

         if( pos == arr.end() || *pos != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 6\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 || arr[5] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the beginning of a small array (x 2 3 4 5 6 7 8)
      {
         test_ = "SmallArray::insert( Iterator, Type&& ) (x 2 3 4 5 6 7 8)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 2, 3, 4, 5, 6, 7, 8 } );

         auto pos = arr.insert( arr.begin(), 1 );

         checkSize    ( arr, 8UL );
         checkCapacity( arr, 8UL );
         checkCount   ( 8U );

         if( pos == arr.end() || *pos != 1 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 1\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ||
             arr[4] != 5 || arr[5] != 6 || arr[6] != 7 || arr[7] != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting in the middle of a small array (1 x 3 4 5 6 7 8 )
      {
         test_ = "SmallArray::insert( Iterator, Type&& ) (1 x 3 4 5 6 7 8)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 3, 4, 5, 6, 7, 8 } );

         auto pos = arr.insert( arr.begin()+1UL, 2 );

         checkSize    ( arr, 8UL );
         checkCapacity( arr, 8UL );
         checkCount   ( 8U );

         if( pos == arr.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ||
             arr[4] != 5 || arr[5] != 6 || arr[6] != 7 || arr[7] != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Inserting at the end of a small array (1 2 3 4 5 6 7 x)
      {
         test_ = "SmallArray::insert( Iterator, Type&& ) (1 2 3 4 5 6 7 x)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5, 6, 7 } );

         auto pos = arr.insert( arr.end(), 8 );

         checkSize    ( arr, 8UL );
         checkCapacity( arr, 8UL );
         checkCount   ( 8U );

         if( pos == arr.end() || *pos != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 8\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ||
             arr[4] != 5 || arr[5] != 6 || arr[6] != 7 || arr[7] != 8 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Inserting an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the SmallArray class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the SmallArray class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testErase()
{
   {
      // Erasing from the beginning of a small array
      {
         test_ = "SmallArray::erase( Iterator ) (x 2 3 4)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4 } );

         auto pos = arr.erase( arr.begin() );

         checkSize    ( arr, 3UL );
         checkCapacity( arr, 3UL );
         checkCount   ( 3U );

         if( pos == arr.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 2 || arr[1] != 3 || arr[2] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the middle of a small array
      {
         test_ = "SmallArray::erase( Iterator ) (1 x 3 4)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4 } );

         auto pos = arr.erase( arr.begin()+1UL );

         checkSize    ( arr, 3UL );
         checkCapacity( arr, 3UL );
         checkCount   ( 3U );

         if( pos == arr.end() || *pos != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 3 || arr[2] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the end of a small array
      {
         test_ = "SmallArray::erase( Iterator ) (1 2 3 x)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4 } );

         auto pos = arr.erase( arr.begin()+3UL );

         checkSize    ( arr, 3UL );
         checkCapacity( arr, 3UL );
         checkCount   ( 3U );

         if( pos != arr.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the beginning of a small array
      {
         test_ = "SmallArray::erase( Iterator ) (x 2 3 4 5 6)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5, 6 } );

         auto pos = arr.erase( arr.begin() );

         checkSize    ( arr, 5UL );
         checkCapacity( arr, 5UL );
         checkCount   ( 5U );

         if( pos == arr.end() || *pos != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 2\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 2 || arr[1] != 3 || arr[2] != 4 || arr[3] != 5 || arr[4] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 2 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the middle of a small array
      {
         test_ = "SmallArray::erase( Iterator ) (1 2 x 4 5 6)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5, 6 } );

         auto pos = arr.erase( arr.begin()+2UL );

         checkSize    ( arr, 5UL );
         checkCapacity( arr, 5UL );
         checkCount   ( 5U );

         if( pos == arr.end() || *pos != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 4 || arr[3] != 5 || arr[4] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the end of a small array
      {
         test_ = "SmallArray::erase( Iterator ) (1 2 3 4 5 x)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5, 6 } );

         auto pos = arr.erase( arr.begin()+5UL );

         checkSize    ( arr, 5UL );
         checkCapacity( arr, 5UL );
         checkCount   ( 5U );

         if( pos != arr.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 || arr[4] != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 5 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      // Erasing from the beginning of a small array
      {
         test_ = "SmallArray::erase( Iterator, Iterator ) (x x 3 4)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4 } );

         auto pos = arr.erase( arr.begin(), arr.begin()+2UL );

         checkSize    ( arr, 2UL );
         checkCapacity( arr, 2UL );
         checkCount   ( 2U );

         if( pos == arr.end() || *pos != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 3 || arr[1] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the middle of a small array
      {
         test_ = "SmallArray::erase( Iterator, Iterator ) (1 x x 4)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4 } );

         auto pos = arr.erase( arr.begin()+1UL, arr.begin()+3UL );

         checkSize    ( arr, 2UL );
         checkCapacity( arr, 2UL );
         checkCount   ( 2U );

         if( pos == arr.end() || *pos != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the end of a small array
      {
         test_ = "SmallArray::erase( Iterator, Iterator ) (1 2 x x)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4 } );

         auto pos = arr.erase( arr.begin()+2UL, arr.begin()+4UL );

         checkSize    ( arr, 2UL );
         checkCapacity( arr, 2UL );
         checkCount   ( 2U );

         if( pos != arr.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the beginning of a small array
      {
         test_ = "SmallArray::erase( Iterator, Iterator ) (x x 3 4 5 6)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5, 6 } );

         auto pos = arr.erase( arr.begin(), arr.begin()+2UL );

         checkSize    ( arr, 4UL );
         checkCapacity( arr, 4UL );
         checkCount   ( 4U );

         if( pos == arr.end() || *pos != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 3\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 3 || arr[1] != 4 || arr[2] != 5 || arr[3] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 3 4 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the middle of a small array
      {
         test_ = "SmallArray::erase( Iterator, Iterator ) (1 2 x x 5 6)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5, 6 } );

         auto pos = arr.erase( arr.begin()+2UL, arr.begin()+4UL );

         checkSize    ( arr, 4UL );
         checkCapacity( arr, 4UL );
         checkCount   ( 4U );

         if( pos == arr.end() || *pos != 5 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << *pos << "\n"
                << "   Expected value: 5\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 5 || arr[3] != 6 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 5 6 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing from the end of a small array
      {
         test_ = "SmallArray::erase( Iterator, Iterator ) (1 2 3 4 x x)";

         blaze::SmallArray<IntResource,N> arr( InitList{ 1, 2, 3, 4, 5, 6 } );

         auto pos = arr.erase( arr.begin()+4UL, arr.begin()+6UL );

         checkSize    ( arr, 4UL );
         checkCapacity( arr, 4UL );
         checkCount   ( 4U );

         if( pos != arr.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( arr[0] != 1 || arr[1] != 2 || arr[2] != 3 || arr[3] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an element failed\n"
                << " Details:\n"
                << "   Result:\n" << arr << "\n"
                << "   Expected result:\n( 1 2 3 4 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the SmallArray class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the SmallArray class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::testSwap()
{
   {
      test_ = "SmallArray swap (size 3 and size 4)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3 } );
      blaze::SmallArray<IntResource,N> arr2( InitList{ 4, 3, 2, 1 } );

      swap( arr1, arr2 );

      checkSize    ( arr1, 4UL );
      checkCapacity( arr1, 4UL );
      checkCount   ( 7U );

      if( arr1[0] != 4 || arr1[1] != 3 || arr1[2] != 2 || arr1[3] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first array failed\n"
             << " Details:\n"
             << "   Result:\n" << arr1 << "\n"
             << "   Expected result:\n( 4 3 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkSize    ( arr2, 3UL );
      checkCapacity( arr2, 3UL );
      checkCount   ( 7U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second array failed\n"
             << " Details:\n"
             << "   Result:\n" << arr1 << "\n"
             << "   Expected result:\n( 1 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray swap (size 3 and size 7)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3 } );
      blaze::SmallArray<IntResource,N> arr2( InitList{ 7, 6, 5, 4, 3, 2, 1 } );

      swap( arr1, arr2 );

      checkSize    ( arr1, 7UL );
      checkCapacity( arr1, 7UL );
      checkCount   ( 10U );

      if( arr1[0] != 7 || arr1[1] != 6 || arr1[2] != 5 || arr1[3] != 4 ||
          arr1[4] != 3 || arr1[5] != 2 || arr1[6] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first array failed\n"
             << " Details:\n"
             << "   Result:\n" << arr1 << "\n"
             << "   Expected result:\n( 7 6 5 4 3 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkSize    ( arr2, 3UL );
      checkCapacity( arr2, 3UL );
      checkCount   ( 10U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second array failed\n"
             << " Details:\n"
             << "   Result:\n" << arr1 << "\n"
             << "   Expected result:\n( 1 2 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray swap (size 8 and size 4)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4, 5, 6, 7, 8 } );
      blaze::SmallArray<IntResource,N> arr2( InitList{ 4, 3, 2, 1 } );

      swap( arr1, arr2 );

      checkSize    ( arr1, 4UL );
      checkCapacity( arr1, 4UL );
      checkCount   ( 12U );

      if( arr1[0] != 4 || arr1[1] != 3 || arr1[2] != 2 || arr1[3] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first array failed\n"
             << " Details:\n"
             << "   Result:\n" << arr1 << "\n"
             << "   Expected result:\n( 4 3 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkSize    ( arr2, 8UL );
      checkCapacity( arr2, 8UL );
      checkCount   ( 12U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 ||
          arr2[4] != 5 || arr2[5] != 6 || arr2[6] != 7 || arr2[7] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second array failed\n"
             << " Details:\n"
             << "   Result:\n" << arr1 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "SmallArray swap (size 8 and size 7)";

      blaze::SmallArray<IntResource,N> arr1( InitList{ 1, 2, 3, 4, 5, 6, 7, 8 } );
      blaze::SmallArray<IntResource,N> arr2( InitList{ 7, 6, 5, 4, 3, 2, 1 } );

      swap( arr1, arr2 );

      checkSize    ( arr1, 7UL );
      checkCapacity( arr1, 7UL );
      checkCount   ( 15U );

      if( arr1[0] != 7 || arr1[1] != 6 || arr1[2] != 5 || arr1[3] != 4 ||
          arr1[4] != 3 || arr1[5] != 2 || arr1[6] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the first array failed\n"
             << " Details:\n"
             << "   Result:\n" << arr1 << "\n"
             << "   Expected result:\n( 7 6 5 4 3 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }

      checkSize    ( arr2, 8UL );
      checkCapacity( arr2, 8UL );
      checkCount   ( 15U );

      if( arr2[0] != 1 || arr2[1] != 2 || arr2[2] != 3 || arr2[3] != 4 ||
          arr2[4] != 5 || arr2[5] != 6 || arr2[6] != 7 || arr2[7] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Swapping the second array failed\n"
             << " Details:\n"
             << "   Result:\n" << arr1 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 6 7 8 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the size of the given small array.
//
// \param array The small array to be checked.
// \param expectedSize The expected size of the small array.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the size of the given small array. In case the actual size does not
// correspond to the given expected size, a \a std::runtime_error exception is thrown.
*/
template< size_t N >       // Number of preallocated elements
template< typename Type >  // Type of the small array
void ClassTest<N>::checkSize( const Type& array, size_t expectedSize ) const
{
   if( array.size() != expectedSize ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid size detected\n"
          << " Details:\n"
          << "   Size         : " << array.size() << "\n"
          << "   Expected size: " << expectedSize << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the capacity of the given small array.
//
// \param array The small array to be checked.
// \param minCapacity The expected minimum capacity of the small array.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of the given small array. In case the actual capacity
// is smaller than the given expected minimum capacity, a \a std::runtime_error exception is
// thrown.
*/
template< size_t N >       // Number of preallocated elements
template< typename Type >  // Type of the small array
void ClassTest<N>::checkCapacity( const Type& array, size_t minCapacity ) const
{
   if( array.capacity() < minCapacity ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Capacity                 : " << array.capacity() << "\n"
          << "   Expected minimum capacity: " << minCapacity << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of 'IntResource' instances.
//
// \param expectedCount The expected number of 'IntResource' instances.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the current number of instances of the 'IntResource' class. In case the
// actual number does not correspond to the actual number, a \a std::runtime_error exception is
// thrown.
*/
template< size_t N >  // Number of preallocated elements
void ClassTest<N>::checkCount( size_t expectedCount ) const
{
   if( IntResource::getCount() != expectedCount ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid count detected\n"
          << " Details:\n"
          << "   Count         : " << IntResource::getCount() << "\n"
          << "   Expected count: " << expectedCount << "\n";
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
/*!\brief Testing the functionality of the SmallArray class template.
//
// \return void
*/
inline void runTest()
{
   ClassTest<0UL>();
   ClassTest<4UL>();
   ClassTest<5UL>();
   ClassTest<6UL>();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the SmallArray class test.
*/
#define RUN_SMALLARRAY_CLASS_TEST \
   blazetest::utiltest::smallarray::runTest();
/*! \endcond */
//*************************************************************************************************

} // namespace smallarray

} // namespace utiltest

} // namespace blazetest

#endif
