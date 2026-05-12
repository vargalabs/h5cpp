//=================================================================================================
/*!
//  \file blazetest/utiltest/alignedallocator/ClassTest.h
//  \brief Header file for the AlignedAllocator test
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

#ifndef _BLAZETEST_UTILTEST_ALIGNEDALLOCATOR_CLASSTEST_H_
#define _BLAZETEST_UTILTEST_ALIGNEDALLOCATOR_CLASSTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <typeinfo>
#include <blaze/util/AlignedAllocator.h>
#include <blaze/util/AlignedArray.h>
#include <blaze/util/AlignmentCheck.h>
#include <blaze/util/typetraits/AlignmentOf.h>


namespace blazetest {

namespace utiltest {

namespace alignedallocator {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for the test of the AlignedAllocator class template.
//
// This class represents the collection of tests for the AlignedAllocator class template.
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
   //**Private class Aligned16*********************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief 16-byte aligned helper class.
   */
   struct Aligned16
   {
      blaze::AlignedArray<int,16UL,16UL> array_;
   };
   /*! \endcond */
   //**********************************************************************************************

   //**Private class Aligned16*********************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief 32-byte aligned helper class.
   */
   struct Aligned32
   {
      blaze::AlignedArray<int,16UL,32UL> array_;
   };
   /*! \endcond */
   //**********************************************************************************************

   //**Private class Aligned16*********************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief 64-byte aligned helper class.
   */
   struct Aligned64
   {
      blaze::AlignedArray<int,16UL,64UL> array_;
   };
   /*! \endcond */
   //**********************************************************************************************

   //**Private class Aligned128********************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief 128-byte aligned helper class.
   */
   struct Aligned128
   {
      blaze::AlignedArray<int,16UL,128UL> array_;
   };
   /*! \endcond */
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   template< typename T > size_t getAlignment( T* ptr ) const;
   //@}
   //**********************************************************************************************

   //**Test functions******************************************************************************
   /*!\name Test functions */
   //@{
   template< typename T > void testAllocation();
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
/*!\brief Test of the allocation/deallocation of for a specific data type.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs an allocation/deallocation of aligned memory for the given type \a T.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >
void ClassTest::testAllocation()
{
   const size_t numObjects( 64UL );

   blaze::AlignedAllocator<T> allocator;
   T* const ptr = allocator.allocate( numObjects );

   if( !blaze::checkAlignment( ptr ) ) {
      std::ostringstream oss;
      oss << " Test: Allocation test for type '" << typeid( T ).name() << "'\n"
          << " Error: Invalid alignment detected\n"
          << " Details:\n"
          << "   Detected alignment = " << getAlignment( ptr ) << "-bit\n"
          << "   Expected alignment = " << blaze::AlignmentOf<T>::value << "-bit\n";
      throw std::runtime_error( oss.str() );
   }

   allocator.deallocate( ptr, numObjects );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Determine the alignment of the given address.
//
// \param address The given address.
// \return The alignment of the given address.
*/
template< typename T >
size_t ClassTest::getAlignment( T* address ) const
{
   size_t alignment( 2UL );

   for( ; alignment<2048UL; alignment*=2UL )
   {
      if( !( reinterpret_cast<size_t>( address ) % alignment ) ) {
         alignment /= 2UL;
         break;
      }
   }

   return alignment;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the functionality of the AlignedAllocator class template.
//
// \return void
*/
inline void runTest()
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
/*!\brief Macro for the execution of the AlignedAllocator class test.
*/
#define RUN_ALIGNEDALLOCATOR_CLASS_TEST \
   blazetest::utiltest::alignedallocator::runTest();
/*! \endcond */
//*************************************************************************************************

} // namespace alignedallocator

} // namespace utiltest

} // namespace blazetest

#endif
