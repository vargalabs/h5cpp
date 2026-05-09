//=================================================================================================
/*!
//  \file src/main/CG.cpp
//  \brief Source file for the conjugate gradient benchmark
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

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/Random.h>
#include <blaze/util/Timing.h>
#include <blazemark/blaze/CG.h>
#include <blazemark/boost/CG.h>
#include <blazemark/eigen/CG.h>
#include <blazemark/gmm/CG.h>
#include <blazemark/mtl/CG.h>
#include <blazemark/system/Boost.h>
#include <blazemark/system/Config.h>
#include <blazemark/system/Eigen.h>
#include <blazemark/system/GMM.h>
#include <blazemark/system/MTL.h>
#include <blazemark/system/Types.h>
#include <blazemark/util/Benchmarks.h>
#include <blazemark/util/Parser.h>
#include <blazemark/util/SolverRun.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


//*************************************************************************************************
// Using declarations
//*************************************************************************************************

using blazemark::Benchmarks;
using blazemark::Parser;
using blazemark::SolverRun;




//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Type of a benchmark run.
//
// This type definition specifies the type of a single benchmark run for the conjugate gradient
// benchmark.
*/
using Run = SolverRun;
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Estimating the necessary number of steps and iterations for each benchmark.
//
// \param run The parameters for the benchmark run.
// \return void
//
// This function estimates the necessary number of steps for the given benchmark based on the
// performance of the Blaze library.
*/
void estimateSteps( Run& run )
{
   using blazemark::element_t;
   using blaze::columnVector;
   using blaze::rowMajor;

   ::blaze::setSeed( ::blazemark::seed );

   const size_t N ( run.getSize() );
   const size_t NN( N*N );

   size_t iterations( run.getIterations() );
   if( iterations == 0UL || iterations > NN ) {
      iterations = NN;
   }

   std::vector<size_t> nnz( NN, 5UL );
   for( size_t i=0UL; i<N; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         if( i == 0UL || i == N-1UL ) --nnz[i*N+j];
         if( j == 0UL || j == N-1UL ) --nnz[i*N+j];
      }
   }

   ::blaze::CompressedMatrix<element_t,rowMajor> A( NN, NN, nnz );
   ::blaze::DynamicVector<element_t,columnVector> x( NN ), b( NN, 0 ), r( NN ), d( NN ), h( NN );
   element_t alpha, beta, delta;
   size_t iteration( 0UL );
   ::blaze::timing::WcTimer timer;

   for( size_t i=0UL; i<N; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         if( i > 0UL   ) A.append( i*N+j, (i-1UL)*N+j, -1.0 );  // Top neighbor
         if( j > 0UL   ) A.append( i*N+j, i*N+j-1UL  , -1.0 );  // Left neighbor
         A.append( i*N+j, i*N+j, 4.0 );
         if( j < N-1UL ) A.append( i*N+j, i*N+j+1UL  , -1.0 );  // Right neighbor
         if( i < N-1UL ) A.append( i*N+j, (i+1UL)*N+j, -1.0 );  // Bottom neighbor
      }
   }

   ::blaze::setSeed( ::blazemark::seed );

   for( size_t i=0UL; i<NN; ++i ) {
      x[i] = ::blaze::rand<element_t>();
   }

   timer.start();

   r = A * x + b;
   delta = trans(r) * r;
   d = -r;

   for( ; iteration<iterations; ++iteration ) {
      h = A * d;
      alpha = delta / ( trans(d) * h );
      x += alpha * d;
      r += alpha * h;
      beta = trans(r) * r;
      if( std::sqrt( beta ) < 1E-8 ) break;
      d = ( beta / delta ) * d - r;
      delta = beta;
   }

   timer.end();

   if( x.size() != NN )
      std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

   if( timer.last() > blazemark::runtime ) {
      const size_t maxIterations( iteration * ( blazemark::runtime / timer.last() ) );
      iteration = blaze::max( 1UL, maxIterations );
   }
   run.setIterations( iteration );

   if( run.getSteps() == 0UL ) {
      if( timer.last() != 0.0 ) {
         const size_t estimatedSteps( blazemark::runtime / timer.last() );
         run.setSteps( blaze::max( 1UL, estimatedSteps ) );
      }
      else {
         run.setSteps( static_cast<size_t>( blazemark::runtime / 1E-8 ) );
      }
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  BENCHMARK FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conjugate gradient benchmark function.
//
// \param runs The specified benchmark runs.
// \param benchmarks The selection of benchmarks.
// \return void
*/
void cg( std::vector<Run>& runs, Benchmarks benchmarks )
{
   std::cout << std::left;

   std::sort( runs.begin(), runs.end() );

   for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
      estimateSteps( *run );
   }

   if( benchmarks.runBlaze ) {
      std::cout << "   Blaze [MFlop/s]:\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N         ( run->getSize()  );
         const size_t steps     ( run->getSteps() );
         const size_t iterations( run->getIterations() );
         run->setBlazeResult( blazemark::blaze::cg( N, steps, iterations ) );
         const double mflops( ( ( 13UL*N*N - 8UL*N - 1UL ) * steps +
                                ( 19UL*N*N - 8UL*N ) * steps * iterations ) / run->getBlazeResult() / 1E6 );
         std::cout << "     " << std::setw(12) << N << mflops << std::endl;
      }
   }

#if BLAZEMARK_BOOST_MODE
   if( benchmarks.runBoost ) {
      std::cout << "   Boost uBLAS [MFlop/s]:\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N         ( run->getSize()  );
         const size_t steps     ( run->getSteps() );
         const size_t iterations( run->getIterations() );
         run->setBoostResult( blazemark::boost::cg( N, steps, iterations ) );
         const double mflops( ( ( 13UL*N*N - 8UL*N - 1UL ) * steps +
                                ( 19UL*N*N - 8UL*N ) * steps * iterations ) / run->getBoostResult() / 1E6 );
         std::cout << "     " << std::setw(12) << N << mflops << std::endl;
      }
   }
#endif

#if BLAZEMARK_GMM_MODE
   if( benchmarks.runGMM ) {
      std::cout << "   GMM++ [MFlop/s]:\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N         ( run->getSize()  );
         const size_t steps     ( run->getSteps() );
         const size_t iterations( run->getIterations() );
         run->setGMMResult( blazemark::gmm::cg( N, steps, iterations ) );
         const double mflops( ( ( 13UL*N*N - 8UL*N - 1UL ) * steps +
                                ( 19UL*N*N - 8UL*N ) * steps * iterations ) / run->getGMMResult() / 1E6 );
         std::cout << "     " << std::setw(12) << N << mflops << std::endl;
      }
   }
#endif

#if BLAZEMARK_MTL_MODE
   if( benchmarks.runMTL ) {
      std::cout << "   MTL [MFlop/s]:\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N         ( run->getSize()  );
         const size_t steps     ( run->getSteps() );
         const size_t iterations( run->getIterations() );
         run->setMTLResult( blazemark::mtl::cg( N, steps, iterations ) );
         const double mflops( ( ( 13UL*N*N - 8UL*N - 1UL ) * steps +
                                ( 19UL*N*N - 8UL*N ) * steps * iterations ) / run->getMTLResult() / 1E6 );
         std::cout << "     " << std::setw(12) << N << mflops << std::endl;
      }
   }
#endif

#if BLAZEMARK_EIGEN_MODE
   if( benchmarks.runEigen ) {
      std::cout << "   Eigen [MFlop/s]:\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N         ( run->getSize()  );
         const size_t steps     ( run->getSteps() );
         const size_t iterations( run->getIterations() );
         run->setEigenResult( blazemark::eigen::cg( N, steps, iterations ) );
         const double mflops( ( ( 13UL*N*N - 8UL*N - 1UL ) * steps +
                                ( 19UL*N*N - 8UL*N ) * steps * iterations ) / run->getEigenResult() / 1E6 );
         std::cout << "     " << std::setw(12) << N << mflops << std::endl;
      }
   }
#endif

   for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
      std::cout << *run;
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The main function for the conjugate gradient benchmark.
//
// \param argc The total number of command line arguments.
// \param argv The array of command line arguments.
// \return void
*/
int main( int argc, char** argv )
{
   std::cout << "\n Conjugate Gradient Method:\n";

   Benchmarks benchmarks;

   try {
      parseCommandLineArguments( argc, argv, benchmarks );
   }
   catch( std::exception& ex ) {
      std::cerr << "   " << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   const std::string installPath( INSTALL_PATH );
   const std::string parameterFile( installPath + "/params/cg.prm" );
   Parser<Run> parser;
   std::vector<Run> runs;

   try {
      parser.parse( parameterFile.c_str(), runs );
   }
   catch( std::exception& ex ) {
      std::cerr << "   Error during parameter extraction: " << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   try {
      cg( runs, benchmarks );
   }
   catch( std::exception& ex ) {
      std::cerr << "   Error during benchmark execution: " << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
