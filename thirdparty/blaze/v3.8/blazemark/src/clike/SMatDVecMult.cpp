//=================================================================================================
/*!
//  \file src/clike/SMatDVecMult.cpp
//  \brief Source file for the C-like sparse matrix/dense vector multiplication kernel
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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iostream>
#include <blaze/util/Indices.h>
#include <blaze/util/Random.h>
#include <blaze/util/Timing.h>
#include <blazemark/clike/SMatDVecMult.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace clike {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief C-like sparse matrix/dense vector multiplication kernel.
//
// \param N The number of rows and columns of the matrix and the size of the vector.
// \param F The number of non-zero elements in each row of the sparse matrix.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the sparse matrix/dense vector multiplication by means of
// a C-like implementation.
*/
double smatdvecmult( size_t N, size_t F, size_t steps )
{
   using ::blazemark::element_t;

   ::blaze::setSeed( seed );

   element_t* value = new element_t[F*N];
   size_t* index = new size_t[F*N];
   size_t* row = new size_t[N+1UL];
   element_t* a = new element_t[N];
   element_t* b = new element_t[N];
   ::blaze::timing::WcTimer timer;
   size_t counter( 0 );

   for( size_t i=0UL; i<N; ++i ) {
      ::blaze::Indices indices( 0UL, N-1UL, F );
      for( ::blaze::Indices::ConstIterator it=indices.begin(); it!=indices.end(); ++it ) {
         value[counter] = ::blaze::rand<element_t>();
         index[counter] = *it;
         ++counter;
      }
      row[i] = i*F;
   }
   row[N] = N*F;

   for( size_t i=0UL; i<N; ++i ) {
      a[i] = ::blaze::rand<element_t>();
   }

   for( size_t i=0UL; i<N; ++i ) {
      b[i] = element_t(0);
      const size_t begin( row[i    ] );
      const size_t end  ( row[i+1UL] );
      for( size_t j=begin; j!=end; ++j ) {
         b[i] += value[j] * a[index[j]];
      }
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         for( size_t i=0UL; i<N; ++i ) {
            b[i] = element_t(0);
            const size_t begin( row[i    ] );
            const size_t end  ( row[i+1UL] );
            for( size_t j=begin; j!=end; ++j ) {
               b[i] += value[j] * a[index[j]];
            }
         }
      }
      timer.end();

      if( b[0] < element_t(0) )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   delete[] value;
   delete[] index;
   delete[] row;
   delete[] a;
   delete[] b;

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " C-like kernel 'smatdvecmult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace clike

} // namespace blazemark
