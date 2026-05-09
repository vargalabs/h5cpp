//=================================================================================================
/*!
//  \file src/mtl/TMat6Vec6Mult.cpp
//  \brief Source file for the MTL 6D transpose matrix/vector multiplication kernel
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
#include <boost/numeric/mtl/mtl.hpp>
#include <blaze/util/Timing.h>
#include <blazemark/mtl/init/Dense2D.h>
#include <blazemark/mtl/init/DenseVector.h>
#include <blazemark/mtl/TMat6Vec6Mult.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace mtl {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief MTL 6-dimensional transpose matrix/vector multiplication kernel.
//
// \param N The number of 6D vectors to be computed.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the 6-dimensional transpose matrix/vector multiplication by
// means of the MTL functionality.
*/
double tmat6vec6mult( size_t N, size_t steps )
{
   using ::blazemark::element_t;

   using col_major    = ::mtl::tag::col_major;
   using parameters   = ::mtl::mat::parameters<col_major>;
   using dense2D      = ::mtl::dense2D<element_t,parameters>;
   using dense_vector = ::mtl::dense_vector<element_t>;

   ::blaze::setSeed( seed );

   ::std::vector<dense2D> A( N );
   ::std::vector<dense_vector> a( N ), b( N );
   ::blaze::timing::WcTimer timer;

   for( size_t i=0UL; i<N; ++i ) {
      A[i].change_dim( 6UL, 6UL );
      a[i].change_dim( 6UL );
      b[i].change_dim( 6UL );
      init( A[i] );
      init( a[i] );
   }

   for( size_t i=0UL; i<N; ++i ) {
      b[i] = A[i] * a[i];
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL, i=0UL; step<steps; ++step, ++i ) {
         if( i == N ) i = 0UL;
         b[i] = A[i] * a[i];
      }
      timer.end();

      for( size_t i=0UL; i<N; ++i )
         if( b[i][0] < element_t(0) )
            std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " MTL kernel 'tmat6vec6mult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace mtl

} // namespace blazemark
