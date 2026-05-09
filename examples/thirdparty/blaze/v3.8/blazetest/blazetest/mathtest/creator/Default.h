//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/Default.h
//  \brief Header file for the Creator class template
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_DEFAULT_H_
#define _BLAZETEST_MATHTEST_CREATOR_DEFAULT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/Boolean.h>
#include <blaze/util/constraints/Builtin.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Void.h>
#include <blaze/util/constraints/Volatile.h>
#include <blazetest/mathtest/creator/Policies.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default creator for random built-in data values.
//
// The Creator class creates random values of the given data type \a T based on the given
// creation policy \a P.
*/
template< typename T >  // Type to be created
class Creator
{
 public:
   //**Type definitions****************************************************************************
   using Type = T;  //!< Type to be created by the Creator.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   // No explicitly declared constructor.
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   // No explicitly declared copy assignment operator.
   inline T operator()() const;
   template< typename CP > inline T operator()( const CP& policy ) const;
   //@}
   //**********************************************************************************************

 private:
   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE    ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST       ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE    ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOID        ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_BOOLEAN_TYPE( T );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created built-in value.
//
// \return The randomly generated built-in value.
*/
template< typename T >  // Type to be created
inline T Creator<T>::operator()() const
{
   return (*this)( Default() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a randomly created built-in value.
//
// \param policy The creation policy for the built-in value.
// \return The randomly generated built-in value.
*/
template< typename T >   // Type to be created
template< typename CP >  // Creation policy
inline T Creator<T>::operator()( const CP& policy ) const
{
   return policy.template create<T>();
}
//*************************************************************************************************

} // namespace blazetest

#endif
