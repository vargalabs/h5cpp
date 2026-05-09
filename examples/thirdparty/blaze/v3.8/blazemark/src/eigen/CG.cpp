//=================================================================================================
/*!
//  \file src/eigen/CG.cpp
//  \brief Source file for the Eigen conjugate gradient kernel
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
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <blaze/util/NumericCast.h>
#include <blaze/util/Timing.h>
#include <blazemark/eigen/CG.h>
#include <blazemark/eigen/init/Matrix.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace eigen {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Eigen conjugate gradient kernel.
//
// \param N The number of rows and columns of the 2D discretized grid.
// \param steps The number of solving steps to perform.
// \param iterations The number of iterations to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the conjugate gradient method by means of the Eigen
// functionality.
*/
double cg( size_t N, size_t steps, size_t iterations )
{
   using ::blazemark::element_t;
   using ::blaze::numeric_cast;
   using ::Eigen::Dynamic;
   using ::Eigen::RowMajor;

   ::blaze::setSeed( seed );

   const size_t NN( N*N );

   ::Eigen::SparseMatrix<element_t,RowMajor,EigenSparseIndexType> A( NN, NN );
   ::Eigen::Matrix<element_t,Dynamic,1> x( NN ), b( NN ), r( NN ), d( NN ), h( NN ), start( NN );
   element_t alpha, beta, delta;
   ::blaze::timing::WcTimer timer;

   A.reserve( NN*5UL );

   for( size_t i=0UL; i<N; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         A.startVec( i*N+j );
         if( i > 0UL   ) A.insertBack(i*N+j,(i-1UL)*N+j) = -1.0;  // Top neighbor
         if( j > 0UL   ) A.insertBack(i*N+j,i*N+j-1UL  ) = -1.0;  // Left neighbor
         A.insertBack(i*N+j,i*N+j) = 4.0;
         if( j < N-1UL ) A.insertBack(i*N+j,i*N+j+1UL  ) = -1.0;  // Right neighbor
         if( i < N-1UL ) A.insertBack(i*N+j,(i+1UL)*N+j) = -1.0;  // Bottom neighbor
      }
   }

   A.finalize();

   for( size_t i=0UL; i<NN; ++i ) {
      b[i] = element_t(0);
   }

   init( start );

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step )
      {
         x.noalias() = start;
         r.noalias() = A * x - b;
         delta = r.transpose() * r;
         d.noalias() = -r;

         for( size_t iteration=0UL; iteration<iterations; ++iteration )
         {
            h = A * d;
            alpha = delta / ( d.transpose() * h );
            x += alpha * d;
            r += alpha * h;
            beta = r.transpose() * r;
            d.noalias() = ( beta / delta ) * d - r;
            delta = beta;
         }
      }

      timer.end();

      if( numeric_cast<size_t>( x.size() ) != NN )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Eigen kernel 'cg': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace eigen

} // namespace blazemark
