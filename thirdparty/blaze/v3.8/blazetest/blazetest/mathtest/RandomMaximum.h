//=================================================================================================
/*!
//  \file blazetest/mathtest/RandomMaximum.h
//  \brief Header file for the RandomMaximum class template
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

#ifndef _BLAZETEST_MATHTEST_RANDOMMAXIMUM_H_
#define _BLAZETEST_MATHTEST_RANDOMMAXIMUM_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <new>
#include <blazetest/mathtest/RandomLimits.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS TEMPLATE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Maximum value for random initializations.
//
// The RandomMaximum class is a wrapper class around the functionality of the RandomLimits class
// template. It represents the largest allowed value for a random initialization of a numeric
// data type. The RandomMaximum class can be implicitly converted to all built-in, numeric data
// types.\n
// In order to handle maximum random values conveniently, the global RandomMaximum instance
// blazetest::randmax is provided, which can be used wherever a numeric built-in data type is
// required:

   \code
   int    i = randmax;  // Assigns the largest allowed signed integer value
   float  f = randmax;  // Assigns the largest allowed single precision value
   double d = randmax;  // Assigns the largest allowed double precision value
   \endcode
*/
class RandomMaximum
{
 public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline RandomMaximum();
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   template< typename T >
   inline operator const T() const;
   //@}
   //**********************************************************************************************

 private:
   //**Forbidden operations************************************************************************
   /*!\name Forbidden operations */
   //@{
   RandomMaximum& operator=( const RandomMaximum& );
   void* operator&() const;

   void* operator new  ( std::size_t ) /*throw( std::bad_alloc )*/;
   void* operator new[]( std::size_t ) /*throw( std::bad_alloc )*/;
   void* operator new  ( std::size_t, const std::nothrow_t& ) /*throw()*/;
   void* operator new[]( std::size_t, const std::nothrow_t& ) /*throw()*/;

   void operator delete  ( void* ) /*throw()*/;
   void operator delete[]( void* ) /*throw()*/;
   void operator delete  ( void*, const std::nothrow_t& ) /*throw()*/;
   void operator delete[]( void*, const std::nothrow_t& ) /*throw()*/;
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
/*!\brief The default constructor of the RandomMaximum class.
*/
inline RandomMaximum::RandomMaximum()
{}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion operator to the required data type.
//
// The conversion operator returns the largest allowed random value for the given data type \a T.
*/
template< typename T >
inline RandomMaximum::operator const T() const
{
   return RandomLimits<T>::max();
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RANDOMMAXIMUM VALUE
//
//=================================================================================================

//*************************************************************************************************
const RandomMaximum randmax;
//*************************************************************************************************

} // namespace blazetest

#endif
