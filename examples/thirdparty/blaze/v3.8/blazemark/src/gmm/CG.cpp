//=================================================================================================
/*!
//  \file src/gmm/CG.cpp
//  \brief Source file for the GMM++ conjugate gradient kernel
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
#include <gmm/gmm.h>
#include <blaze/util/Timing.h>
#include <blazemark/gmm/CG.h>
#include <blazemark/gmm/init/Vector.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace gmm {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief GMM++ conjugate gradient kernel.
//
// \param N The number of rows and columns of the 2D discretized grid.
// \param steps The number of solving steps to perform.
// \param iterations The number of iterations to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the conjugate gradient method by means of the GMM++
// functionality.
*/
double cg( size_t N, size_t steps, size_t iterations )
{
   using ::blazemark::element_t;

   ::blaze::setSeed( seed );

   const size_t NN( N*N );

   ::gmm::row_matrix< ::gmm::wsvector<element_t> > T( NN, NN );
   ::gmm::csr_matrix<element_t> A( NN, NN );
   ::std::vector<element_t> x( NN ), b( NN ), r( NN ), d( NN ), h( NN ), start( NN );
   element_t alpha, beta, delta;
   ::blaze::timing::WcTimer timer;

   for( size_t i=0UL; i<N; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         if( i > 0UL   ) T(i*N+j,(i-1UL)*N+j) = -1.0;  // Top neighbor
         if( j > 0UL   ) T(i*N+j,i*N+j-1UL  ) = -1.0;  // Left neighbor
         T(i*N+j,i*N+j) = 4.0;
         if( j < N-1UL ) T(i*N+j,i*N+j+1UL  ) = -1.0;  // Right neighbor
         if( i < N-1UL ) T(i*N+j,(i+1UL)*N+j) = -1.0;  // Bottom neighbor
      }
   }

   copy( T, A );

   ::gmm::clear( b );
   init( start );

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step )
      {
         x = start;
         mult( A, x, r );
         ::gmm::add( ::gmm::scaled( b, element_t(-1) ), r );
         delta = ::gmm::vect_sp( r, r );
         copy( ::gmm::scaled( r, element_t(-1) ), d );

         for( size_t iteration=0UL; iteration<iterations; ++iteration )
         {
            mult( A, d, h );
            alpha = delta / ::gmm::vect_sp( d, h );
            ::gmm::add( ::gmm::scaled( d, alpha ), x, x );
            ::gmm::add( ::gmm::scaled( h, alpha ), r, r );
            beta = ::gmm::vect_sp( r, r );
            ::gmm::add( ::gmm::scaled( r, element_t(-1) ), ::gmm::scaled( d, beta / delta ), d );
            delta = beta;
         }
      }
      timer.end();

      if( x.size() != NN )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " GMM++ kernel 'cg': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace gmm

} // namespace blazemark
