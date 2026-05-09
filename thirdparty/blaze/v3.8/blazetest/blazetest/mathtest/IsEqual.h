//=================================================================================================
/*!
//  \file blazetest/mathtest/IsEqual.h
//  \brief Header file for the isEqual comparison
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

#ifndef _BLAZETEST_MATHTEST_ISEQUAL_H_
#define _BLAZETEST_MATHTEST_ISEQUAL_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/shims/Equal.h>
#include <blaze/util/typetraits/IsNumeric.h>


namespace blazetest {

//=================================================================================================
//
//  COMPARISON FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary comparison function.
//
// \param m1 The left-hand side element to be compared.
// \param m2 The right-hand side element to be compared.
// \return \a true in case the two elements are equal, \a false in case the elements are not equal.
//
// This function performs an excessive comparison of the two given elements. It utilizes both
// equality and inequality comparison and additionally swaps the two elements.
*/
template< typename T1    // Type of the first element
        , typename T2 >  // Type of the second element
bool isEqual( const T1& m1, const T2& m2 )
{
   if( blaze::IsNumeric<T1>::value && blaze::IsNumeric<T2>::value )
      return blaze::equal( m1, m2 );
   else
      return ( m1 == m2 && !( m1 != m2 ) && m2 == m1 && !( m2 != m1 ) );
}
//*************************************************************************************************

} // namespace blazetest

#endif
