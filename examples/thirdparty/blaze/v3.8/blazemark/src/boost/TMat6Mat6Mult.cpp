//=================================================================================================
/*!
//  \file src/boost/TMat6Mat6Mult.cpp
//  \brief Source file for the Boost 6D transpose matrix/matrix multiplication kernel
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
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <blaze/util/Timing.h>
#include <blazemark/boost/init/Matrix.h>
#include <blazemark/boost/TMat6Mat6Mult.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace boost {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Boost uBLAS 6-dimensional transpose matrix/matrix multiplication kernel.
//
// \param N The number of 6x6 matrices to be computed.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the 6-dimensional transpose matrix/matrix multiplication by
// means of the Boost uBLAS functionality.
*/
double tmat6mat6mult( size_t N, size_t steps )
{
   using ::blazemark::element_t;
   using ::boost::numeric::ublas::row_major;
   using ::boost::numeric::ublas::column_major;

   ::blaze::setSeed( seed );

   ::std::vector< ::boost::numeric::ublas::matrix<element_t,column_major> > A( N ), C( N );
   ::std::vector< ::boost::numeric::ublas::matrix<element_t,row_major> > B( N );
   ::blaze::timing::WcTimer timer;

   for( size_t i=0UL; i<N; ++i ) {
      A[i].resize( 6UL, 6UL );
      B[i].resize( 6UL, 6UL );
      C[i].resize( 6UL, 6UL );
      init( A[i] );
      init( B[i] );
   }

   for( size_t i=0UL; i<N; ++i ) {
      noalias( C[i] ) = prod( A[i], B[i] );
   }

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL, i=0UL; step<steps; ++step, ++i ) {
         if( i == N ) i = 0UL;
         noalias( C[i] ) = prod( A[i], B[i] );
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
      std::cerr << " Boost uBLAS kernel 'tmat6mat6mult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace boost

} // namespace blazemark
