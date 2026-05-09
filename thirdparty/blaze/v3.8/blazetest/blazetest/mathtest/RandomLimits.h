//=================================================================================================
/*!
//  \file blazetest/mathtest/RandomLimits.h
//  \brief Header file for the RandomLimits class template
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

#ifndef _BLAZETEST_MATHTEST_RANDOMLIMITS_H_
#define _BLAZETEST_MATHTEST_RANDOMLIMITS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsSigned.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization limits for all built-in data types
//
// The RandomLimits class provides minimum and maximum limits for the random initialization of
// all numerical built-in data types. Via the \a min and \a max functions all numerical data
// types are restricted to a fixed range of values. Unsigned integral data types are restricted
// to the range \f$ [0..10] \f$, signed integral data types to the range \f$ [-10..10] \f$, and
// floating point data types to the range \f$ [-1..1) \f$.
//
// Code examples:

   \code
   // Smallest allowed signed integer value
   int i = RandomLimits<int>::min();

   // Largest allowed double precision floating point value
   double d = RandomLimits<double>::max();
   \endcode
*/
template< typename Type >
struct RandomLimits
{
   //**Limit functions*****************************************************************************
   /*!\name Limit functions */
   //@{
   static inline Type min() noexcept;
   static inline Type max() noexcept;
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( Type );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  LIMIT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization minimum.
//
// \return The smallest allowed initialization value.
*/
template< typename Type >
inline Type RandomLimits<Type>::min() noexcept
{
   if( blaze::IsFloatingPoint<Type>::value )
      return Type(-1);
   else if( blaze::IsSigned<Type>::value )
      return Type(-10);
   else
      return Type(0);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization maximum.
//
// \return The largest allowed initialization value.
*/
template< typename Type >
inline Type RandomLimits<Type>::max() noexcept
{
   if( blaze::IsFloatingPoint<Type>::value )
      return Type(1);
   else
      return Type(10);
}
//*************************************************************************************************

} // namespace blazetest

#endif
