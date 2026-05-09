//=================================================================================================
/*!
//  \file src/gmm/TMat3Vec3Mult.cpp
//  \brief Source file for the GMM++ 3D transpose matrix/vector multiplication kernel
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
#include <gmm/gmm.h>
#include <blaze/util/Timing.h>
#include <blazemark/gmm/init/DenseMatrix.h>
#include <blazemark/gmm/init/Vector.h>
#include <blazemark/gmm/TDMatDVecMult.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace gmm {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief GMM++ 3-dimensional transpose matrix/vector multiplication kernel.
//
// \param N The number of 3D vectors to be computed.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the 3-dimensional transpose matrix/vector multiplication by
// means of the GMM++ functionality.
*/
double tmat3vec3mult( size_t N, size_t steps )
{
   using ::blazemark::element_t;

   ::blaze::setSeed( seed );

   ::std::vector< ::gmm::dense_matrix<element_t> > A( N );
   ::std::vector< ::std::vector<element_t> > a( N ), b( N );
   ::blaze::timing::WcTimer timer;

   for( size_t i=0UL; i<N; ++i ) {
      ::gmm::resize( A[i], 3UL, 3UL );
      ::gmm::resize( a[i], 3UL );
      ::gmm::resize( b[i], 3UL );
      init( A[i] );
      init( a[i] );
   }

   for( size_t i=0UL; i<N; ++i ) {
      mult( A[i], a[i], b[i] );
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL, i=0UL; step<steps; ++step, ++i ) {
         if( i == N ) i = 0UL;
         mult( A[i], a[i], b[i] );
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
      std::cerr << " GMM++ kernel 'tmat3vec3mult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace gmm

} // namespace blazemark
