//=================================================================================================
/*!
//  \file src/blitz/DVecDVecOuter.cpp
//  \brief Source file for the Blitz++ dense vector/dense vector outer product kernel
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
#include <blitz/array.h>
#include <blaze/util/NumericCast.h>
#include <blaze/util/Timing.h>
#include <blazemark/blitz/DVecDVecOuter.h>
#include <blazemark/blitz/init/Array.h>
#include <blazemark/system/Config.h>


namespace blazemark {

namespace blitz {

//=================================================================================================
//
//  KERNEL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Blitz++ dense vector/dense vector outer product kernel.
//
// \param N The size of the vectors for the outer product.
// \param steps The number of iteration steps to perform.
// \return Minimum runtime of the kernel function.
//
// This kernel function implements the dense vector/dense vector outer product by means of
// the Blitz++ functionality.
*/
double dvecdvecouter( size_t N, size_t steps )
{
   using ::blazemark::element_t;
   using ::blaze::numeric_cast;

   ::blaze::setSeed( seed );

   ::blitz::Array<element_t,1> a( N ), b( N );
   ::blitz::Array<element_t,2> A( N, N );
   ::blitz::firstIndex i;
   ::blitz::secondIndex j;
   ::blaze::timing::WcTimer timer;

   init( a );
   init( b );

   A = a(i) * b(j);

   for( size_t rep=0UL; rep<reps; ++rep )
   {
      timer.start();
      for( size_t step=0UL; step<steps; ++step ) {
         A = a(i) * b(j);
      }
      timer.end();

      if( numeric_cast<size_t>( A.rows() ) != N )
         std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

      if( timer.last() > maxtime )
         break;
   }

   const double minTime( timer.min()     );
   const double avgTime( timer.average() );

   if( minTime * ( 1.0 + deviation*0.01 ) < avgTime )
      std::cerr << " Blitz++ kernel 'dvecdvecouter': Time deviation too large!!!\n";

   return minTime;
}
//*************************************************************************************************

} // namespace blitz

} // namespace blazemark
