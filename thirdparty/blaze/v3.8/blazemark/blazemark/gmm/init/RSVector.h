//=================================================================================================
/*!
//  \file blazemark/gmm/init/RSVector.h
//  \brief Header file for the GMM++ sparse vector initialization functions
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

#ifndef _BLAZEMARK_GMM_INIT_RSVECTOR_H_
#define _BLAZEMARK_GMM_INIT_RSVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <gmm/gmm_vector.h>
#include <blaze/util/Indices.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace gmm {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name GMM++ initialization functions */
//@{
template< typename Type >
void init( ::gmm::rsvector<Type>& v, size_t nonzeros );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given sparse vector.
//
// \param v The sparse vector to be initialized.
// \param nonzeros The number of non-zero elements.
// \return void
//
// This function initializes the given sparse vector with random values. The vector will
// be filled with \a nonzeros non-zero elements, whose indices will be randomly determined.
*/
template< typename Type >  // Data type of the vector
void init( ::gmm::rsvector<Type>& v, size_t nonzeros )
{
   const size_t N( vect_size( v ) );

   ::blaze::Indices indices( 0UL, N-1UL, nonzeros );
   for( ::blaze::Indices::ConstIterator it=indices.begin(); it!=indices.end(); ++it ) {
      v[*it] = ::blaze::rand<Type>( 0, 10 );
   }
}
//*************************************************************************************************

} // namespace gmm

} // namespace blazemark

#endif
