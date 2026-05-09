//=================================================================================================
/*!
//  \file src/flens/Mat6TMat6Mult.cpp
//  \brief Source file for the FLENS 6D matrix/transpose matrix multiplication kernel
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
#include <vector>
#include <flens/flens.cxx>
#include <blaze/util/Timing.h>
#include <blazemark/flens/init/GeMatrix.h>
#include <blazemark/flens/Mat6TMat6Mult.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace flens {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief FLENS 6-dimensional matrix/transpose matrix multiplication kernel.
//
// \param N The number of 6x6 matrices to be computed.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the 6-dimensional matrix/transpose matrix multiplication by
// means of the FLENS functionality.
*/
double mat6tmat6mult( size_t N, size_t steps )
{
   using ::blazemark::element_t;

   using RowGeMatrix = ::flens::GeMatrix< ::flens::FullStorage<element_t,::flens::RowMajor> >;
   using ColGeMatrix = ::flens::GeMatrix< ::flens::FullStorage<element_t,::flens::ColMajor> >;

   ::blaze::setSeed( seed );

   ::std::vector< RowGeMatrix > A( N ), C( N );
   ::std::vector< ColGeMatrix > B( N );
   ::blaze::timing::WcTimer timer;

   for( size_t i=0UL; i<N; ++i ) {
      A[i].resize( 6, 6 );
      B[i].resize( 6, 6 );
      C[i].resize( 6, 6 );
      init( A[i] );
      init( B[i] );
   }

   for( size_t i=0UL; i<N; ++i ) {
      C[i] = A[i] * B[i];
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL, i=0UL; step<steps; ++step, ++i ) {
         if( i == N ) i = 0UL;
         C[i] = A[i] * B[i];
      }
      timer.end();

      for( size_t i=0UL; i<N; ++i )
         if( C[i](0,0) < element_t(0) )
            std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " FLENS kernel 'mat6tmat6mult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace flens

} // namespace blazemark
