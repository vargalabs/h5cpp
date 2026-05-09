//=================================================================================================
/*!
//  \file blazetest/utiltest/InstanceCounter.h
//  \brief Header file for InstanceCounter class template
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

#ifndef _BLAZETEST_UTILTEST_INSTANCECOUNTER_H_
#define _BLAZETEST_UTILTEST_INSTANCECOUNTER_H_


namespace blazetest {

namespace utiltest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for classes requiring an instance counter.
//
// The InstanceCounter class template provides the functionality to count the number of instances
// created from a particular class. The functionality can be used by either deriving publicly
// (providing public access to the instance counter) or non-publicly (hiding access to the
// instance counter) from the InstanceCounter class. The InstanceCounter class template requires
// the deriving class as template argument (CRTP pattern):

   \code
   class Resource : public InstanceCounter<Resource>
   {
      // ...
   };
   \endcode
*/
template< typename T >
class InstanceCounter
{
 protected:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline InstanceCounter() noexcept;
   inline InstanceCounter( const InstanceCounter& ) noexcept;
   inline InstanceCounter( InstanceCounter&& ) noexcept;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~InstanceCounter();
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   InstanceCounter& operator=( const InstanceCounter& ) = default;
   InstanceCounter& operator=( InstanceCounter&& ) = default;
   //@}
   //**********************************************************************************************

 public:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static inline unsigned int getCount() noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   static unsigned int counter_;  // The instance counter.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor of InstanceCounter.
*/
template< typename T >
inline InstanceCounter<T>::InstanceCounter() noexcept
{
   ++counter_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor of InstanceCounter.
*/
template< typename T >
inline InstanceCounter<T>::InstanceCounter( const InstanceCounter& ) noexcept
{
   ++counter_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The move constructor of InstanceCounter.
*/
template< typename T >
inline InstanceCounter<T>::InstanceCounter( InstanceCounter&& ) noexcept
{
   ++counter_;
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor of InstanceCounter.
*/
template< typename T >
inline InstanceCounter<T>::~InstanceCounter()
{
   --counter_;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current count of instances.
//
// \return The current count of instances.
*/
template< typename T >
inline unsigned int InstanceCounter<T>::getCount() noexcept
{
   return counter_;
}
//*************************************************************************************************




//=================================================================================================
//
//  DEFINITION AND INITIALIZATION OF THE STATIC MEMBER VARIABLES
//
//=================================================================================================

template< typename T >
unsigned int InstanceCounter<T>::counter_( 0U );

} // namespace utiltest

} // namespace blazetest

#endif
