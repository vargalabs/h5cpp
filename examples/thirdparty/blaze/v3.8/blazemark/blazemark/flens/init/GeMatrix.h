//=================================================================================================
/*!
//  \file blazemark/flens/init/GeMatrix.h
//  \brief Header file for the FLENS dense matrix initialization functions
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

#ifndef _BLAZEMARK_FLENS_INIT_GEMATRIX_H_
#define _BLAZEMARK_FLENS_INIT_GEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <flens/matrixtypes/general/impl/gematrix.h>
#include <blaze/util/Random.h>
#include <blazemark/system/Types.h>


namespace blazemark {

namespace flens {

//=================================================================================================
//
//  INITIALIZATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name FLENS initialization functions */
//@{
template< typename Type >
void init( ::flens::GeMatrix< ::flens::FullStorage<Type,::flens::RowMajor> >& m );

template< typename Type >
void init( ::flens::GeMatrix< ::flens::FullStorage<Type,::flens::ColMajor> >& m );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given row-major dense matrix.
//
// \param m The row-major dense matrix to be initialized.
// \return void
//
// This function initializes the given row-major dense matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::flens::GeMatrix< ::flens::FullStorage<Type,::flens::RowMajor> >& m )
{
   using IndexType = typename ::flens::GeMatrix< ::flens::FullStorage<Type,::flens::RowMajor> >::IndexType;

   for( IndexType i=m.firstRow(); i<=m.lastRow(); ++i ) {
      for( IndexType j=m.firstCol(); j<=m.lastCol(); ++j ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Random initialization of the given column-major dense matrix.
//
// \param m The column-major dense matrix to be initialized.
// \return void
//
// This function initializes the given column-major dense matrix with random values.
*/
template< typename Type >  // Data type of the matrix
void init( ::flens::GeMatrix< ::flens::FullStorage<Type,::flens::ColMajor> >& m )
{
   using IndexType = typename ::flens::GeMatrix< ::flens::FullStorage<Type,::flens::ColMajor> >::IndexType;

   for( IndexType j=m.firstCol(); j<=m.lastCol(); ++j ) {
      for( IndexType i=m.firstRow(); i<=m.lastRow(); ++i ) {
         m(i,j) = ::blaze::rand<Type>( 0, 10 );
      }
   }
}
//*************************************************************************************************

} // namespace flens

} // namespace blazemark

#endif
