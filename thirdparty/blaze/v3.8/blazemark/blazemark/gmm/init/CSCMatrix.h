//=================================================================================================
/*!
//  \file blazemark/gmm/init/CSCMatrix.h
//  \brief Header file for the GMM++ column-major sparse matrix initialization functions
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

#ifndef _BLAZEMARK_GMM_INIT_CSCMATRIX_H_
#define _BLAZEMARK_GMM_INIT_CSCMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <gmm/gmm_matrix.h>
#include <blaze/util/Indices.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Config.h>
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
void init( ::gmm::csc_matrix<Type>& m, size_t nonzeros );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major sparse matrix.
//
// \param m The column-major sparse matrix to be initialized.
// \param nonzeros The number of non-zero elements per column.
// \return void
//
// This function initializes the given column-major sparse matrix with random values.
// Each column will be filled with \a nonzeros non-zero elements, whose indices are randomly
// determined.
*/
template< typename Type >  // Data type of the matrix
void init( ::gmm::csc_matrix<Type>& m, size_t nonzeros )
{
   const size_t M( mat_nrows(m) );
   const size_t N( mat_ncols(m) );

   ::gmm::col_matrix< ::gmm::wsvector<Type> > tmp( M, N );

   if( structure == band )
   {
      const size_t drange( nonzeros / 2UL );
      const size_t urange( ( nonzeros % 2UL )?( drange ):( drange-1UL ) );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( j >= urange )?( j-urange ):( 0UL ) );
         const size_t iend  ( ( j+drange+1UL < M )?( j+drange+1UL ):( M ) );

         for( size_t i=ibegin; i<iend; ++i ) {
            tmp(i,j) = ::blaze::rand<Type>( 0, 10 );
         }
      }
   }
   else
   {
      for( size_t j=0UL; j<N; ++j ) {
         ::blaze::Indices indices( 0UL, M-1UL, nonzeros );
         for( ::blaze::Indices::ConstIterator it=indices.begin(); it!=indices.end(); ++it ) {
            tmp(*it,j) = ::blaze::rand<Type>( 0, 10 );
         }
      }
   }

   copy( tmp, m );
}
//*************************************************************************************************

} // namespace gmm

} // namespace blazemark

#endif
