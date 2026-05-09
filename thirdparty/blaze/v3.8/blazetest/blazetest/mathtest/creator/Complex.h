//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/Complex.h
//  \brief Specialization of the Creator class template for complex values
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_COMPLEX_H_
#define _BLAZETEST_MATHTEST_CREATOR_COMPLEX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/Boolean.h>
#include <blaze/util/constraints/Builtin.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Void.h>
#include <blaze/util/constraints/Volatile.h>
#include <blazetest/mathtest/creator/Default.h>
#include <blazetest/system/Types.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the Creator class template for complex data types.
//
// This specialization of the Creator class template is able to create random complex values.
*/
template< typename T >  // Type to be created
class Creator< complex<T> >
{
 public:
   //**Type definitions****************************************************************************
   using Type = complex<T>;  //!< Type to be created by the Creator.
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
   const complex<T> operator()() const;
   template< typename CP > const complex<T> operator()( const CP& policy ) const;
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
/*!\brief Returns a randomly created complex value.
//
// \return The randomly generated complex value.
*/
template< typename T >  // Type to be created
inline const complex<T> Creator< complex<T> >::operator()() const
{
   return (*this)( Default() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a randomly created complex value.
//
// \param policy The creation policy for the real and imaginary part.
// \return The randomly generated complex value.
*/
template< typename T >   // Type to be created
template< typename CP >  // Creation policy
inline const complex<T> Creator< complex<T> >::operator()( const CP& policy ) const
{
   return complex<T>( policy.template create<T>(), policy.template create<T>() );
}
//*************************************************************************************************

} // namespace blazetest

#endif
