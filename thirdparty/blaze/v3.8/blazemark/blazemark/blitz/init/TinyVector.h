//=================================================================================================
/*!
//  \file blazemark/blitz/init/TinyVector.h
//  \brief Header file for the Blitz++ tiny vector initialization functions
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

#ifndef _BLAZEMARK_BLITZ_INIT_TINYVECTOR_H_
#define _BLAZEMARK_BLITZ_INIT_TINYVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blitz/tinyvec2.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace blitz {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Blitz++ initialization functions */
//@{
template< typename Type, int N >
void init( ::blitz::TinyVector<Type,N>& v );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given tiny vector.
//
// \param v The tiny vector to be initialized.
// \return void
//
// This function initializes the given tiny vector with random values.
*/
template< typename Type  // Data type of the vector
        , int N >        // Number of elements
void init( ::blitz::TinyVector<Type,N>& v )
{
   for( int i=0; i<N; ++i ) {
      v(i) = ::blaze::rand<Type>( 0, 10 );
   }
}
//*************************************************************************************************

} // namespace blitz

} // namespace blazemark

#endif
