//=================================================================================================
/*!
//  \file src/flens/TDMatDMatMult.cpp
//  \brief Source file for the FLENS transpose dense matrix/dense matrix multiplication kernel
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
#include <flens/flens.cxx>
#include <blaze/util/NumericCast.h>
#include <blaze/util/Timing.h>
#include <blazemark/flens/init/GeMatrix.h>
#include <blazemark/flens/TDMatDMatMult.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace flens {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief FLENS transpose dense matrix/dense matrix multiplication kernel.
//
// \param N The number of rows and columns of the matrices.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the transpose dense matrix/dense matrix multiplication by
// means of the FLENS functionality.
*/
double tdmatdmatmult( size_t N, size_t steps )
{
   using ::blazemark::element_t;
   using ::blaze::numeric_cast;

   using RowGeMatrix = ::flens::GeMatrix< ::flens::FullStorage<element_t,::flens::RowMajor> >;
   using ColGeMatrix = ::flens::GeMatrix< ::flens::FullStorage<element_t,::flens::ColMajor> >;

   ::blaze::setSeed( seed );

   ColGeMatrix A( N, N ), C( N, N );
   RowGeMatrix B( N, N );
   ::blaze::timing::WcTimer timer;

   init( A );
   init( B );

   C = A * B;

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         C = A * B;
      }
      timer.end();

      if( numeric_cast<size_t>( C.numRows() ) != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " FLENS kernel 'tdmatdmatmult': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace flens

} // namespace blazemark
