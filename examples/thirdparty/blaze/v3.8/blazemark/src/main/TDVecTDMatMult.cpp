//=================================================================================================
/*!
//  \file src/main/TDVecDMatMult.cpp
//  \brief Source file for the transpose dense vector/transpose dense matrix multiplication benchmark
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
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/Infinity.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/Random.h>
#include <blaze/util/Timing.h>
#include <blazemark/armadillo/TDVecTDMatMult.h>
#include <blazemark/blas/TDVecTDMatMult.h>
#include <blazemark/blaze/init/DynamicMatrix.h>
#include <blazemark/blaze/init/DynamicVector.h>
#include <blazemark/blaze/TDVecTDMatMult.h>
#include <blazemark/boost/TDVecTDMatMult.h>
#include <blazemark/eigen/TDVecTDMatMult.h>
#include <blazemark/flens/TDVecTDMatMult.h>
#include <blazemark/system/Armadillo.h>
#include <blazemark/system/BLAS.h>
#include <blazemark/system/Boost.h>
#include <blazemark/system/Config.h>
#include <blazemark/system/Eigen.h>
#include <blazemark/system/FLENS.h>
#include <blazemark/system/Types.h>
#include <blazemark/util/Benchmarks.h>
#include <blazemark/util/DynamicDenseRun.h>
#include <blazemark/util/Parser.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


//*************************************************************************************************
// Using declarations
//*************************************************************************************************

using blazemark::Benchmarks;
using blazemark::DynamicDenseRun;
using blazemark::Parser;




//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Type of a benchmark run.
//
// This type definition specifies the type of a single benchmark run for the transpose dense
// vector/transpose dense matrix multiplication benchmark
*/
using Run = DynamicDenseRun;
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Estimating the necessary number of steps for each benchmark.
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
   using blaze::rowVector;
   using blaze::columnMajor;

   ::blaze::setSeed( ::blazemark::seed );

   const size_t N( run.getSize() );

   blaze::DynamicVector<element_t,rowVector> a( N ), b( N );
   blaze::DynamicMatrix<element_t,columnMajor> A( N, N );
   blaze::timing::WcTimer timer;
   double wct( 0.0 );
   size_t steps( 1UL );

   blazemark::blaze::init( a );
   blazemark::blaze::init( A );

   while( true ) {
      timer.start();
      for( size_t i=0UL; i<steps; ++i ) {
         b = a * A;
      }
      timer.end();
      wct = timer.last();
      if( wct >= 0.2 ) break;
      steps *= 2UL;
   }

   if( b.size() != N )
      std::cerr << " Line " << __LINE__ << ": ERROR detected!!!\n";

   const size_t estimatedSteps( ( blazemark::runtime * steps ) / timer.last() );
   run.setSteps( blaze::max( 1UL, estimatedSteps ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimating the necessary number of floating point operations.
//
// \param run The parameters for the benchmark run.
// \return void
//
// This function estimates the number of floating point operations required for a single
// computation of the (composite) arithmetic operation.
*/
void estimateFlops( Run& run )
{
   const size_t N( run.getSize() );

   run.setFlops( 2UL*N*N - N );
}
//*************************************************************************************************




//=================================================================================================
//
//  BENCHMARK FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Transpose dense vector/transpose dense matrix multiplication benchmark function.
//
// \param runs The specified benchmark runs.
// \param benchmarks The selection of benchmarks.
// \return void
*/
void tdvectdmatmult( std::vector<Run>& runs, Benchmarks benchmarks )
{
   std::cout << std::left;

   std::sort( runs.begin(), runs.end() );

   size_t slowSize( blaze::inf );
   for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run )
   {
      estimateFlops( *run );

      if( run->getSteps() == 0UL ) {
         if( run->getSize() < slowSize ) {
            estimateSteps( *run );
            if( run->getSteps() == 1UL )
               slowSize = run->getSize();
         }
         else run->setSteps( 1UL );
      }
   }

#if BLAZEMARK_BLAS_MODE
   if( benchmarks.runBLAS ) {
      std::cout << "   BLAS implementation [MFlop/s]:\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()  );
         const size_t steps( run->getSteps() );
         run->setBLASResult( blazemark::blas::tdvectdmatmult( N, steps ) );
         const double mflops( run->getFlops() * steps / run->getBLASResult() / 1E6 );
         std::cout << "     " << std::setw(12) << N << mflops << std::endl;
      }
   }
#endif

   if( benchmarks.runBlaze ) {
      std::cout << "   Blaze [MFlop/s]:\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()  );
         const size_t steps( run->getSteps() );
         run->setBlazeResult( blazemark::blaze::tdvectdmatmult( N, steps ) );
         const double mflops( run->getFlops() * steps / run->getBlazeResult() / 1E6 );
         std::cout << "     " << std::setw(12) << N << mflops << std::endl;
      }
   }

#if BLAZEMARK_BOOST_MODE
   if( benchmarks.runBoost ) {
      std::cout << "   Boost uBLAS [MFlop/s]:\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()  );
         const size_t steps( run->getSteps() );
         run->setBoostResult( blazemark::boost::tdvectdmatmult( N, steps ) );
         const double mflops( run->getFlops() * steps / run->getBoostResult() / 1E6 );
         std::cout << "     " << std::setw(12) << N << mflops << std::endl;
      }
   }
#endif

#if BLAZEMARK_ARMADILLO_MODE
   if( benchmarks.runArmadillo ) {
      std::cout << "   Armadillo [MFlop/s]:\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()  );
         const size_t steps( run->getSteps() );
         run->setArmadilloResult( blazemark::armadillo::tdvectdmatmult( N, steps ) );
         const double mflops( run->getFlops() * steps / run->getArmadilloResult() / 1E6 );
         std::cout << "     " << std::setw(12) << N << mflops << std::endl;
      }
   }
#endif

#if BLAZEMARK_FLENS_MODE
   if( benchmarks.runFLENS ) {
      std::cout << "   FLENS [MFlop/s]:\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()  );
         const size_t steps( run->getSteps() );
         run->setFLENSResult( blazemark::flens::tdvectdmatmult( N, steps ) );
         const double mflops( run->getFlops() * steps / run->getFLENSResult() / 1E6 );
         std::cout << "     " << std::setw(12) << N << mflops << std::endl;
      }
   }
#endif

#if BLAZEMARK_EIGEN_MODE
   if( benchmarks.runEigen ) {
      std::cout << "   Eigen [MFlop/s]:\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()  );
         const size_t steps( run->getSteps() );
         run->setEigenResult( blazemark::eigen::tdvectdmatmult( N, steps ) );
         const double mflops( run->getFlops() * steps / run->getEigenResult() / 1E6 );
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
/*!\brief The main function for the transpose dense vector/transpose dense matrix multiplication
//        benchmark.
//
// \param argc The total number of command line arguments.
// \param argv The array of command line arguments.
// \return void
*/
int main( int argc, char** argv )
{
   std::cout << "\n Transpose Dense Vector/Transpose Dense Matrix Multiplication:\n";

   Benchmarks benchmarks;

   try {
      parseCommandLineArguments( argc, argv, benchmarks );
   }
   catch( std::exception& ex ) {
      std::cerr << "   " << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   const std::string installPath( INSTALL_PATH );
   const std::string parameterFile( installPath + "/params/tdvectdmatmult.prm" );
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
      tdvectdmatmult( runs, benchmarks );
   }
   catch( std::exception& ex ) {
      std::cerr << "   Error during benchmark execution: " << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
