//=================================================================================================
/*!
//  \file src/main/Custom.cpp
//  \brief Source file for the benchmark for custom expressions
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
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <blaze/math/Infinity.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/Timing.h>
#include <blazemark/armadillo/Custom.h>
#include <blazemark/blaze/Custom.h>
#include <blazemark/blitz/Custom.h>
#include <blazemark/boost/Custom.h>
#include <blazemark/eigen/Custom.h>
#include <blazemark/gmm/Custom.h>
#include <blazemark/mtl/Custom.h>
#include <blazemark/system/Armadillo.h>
#include <blazemark/system/Blitz.h>
#include <blazemark/system/Boost.h>
#include <blazemark/system/Config.h>
#include <blazemark/system/Eigen.h>
#include <blazemark/system/GMM.h>
#include <blazemark/system/MTL.h>
#include <blazemark/system/Types.h>
#include <blazemark/util/Benchmarks.h>
#include <blazemark/util/DynamicSparseRun.h>
#include <blazemark/util/Parser.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


//*************************************************************************************************
// Using declarations
//*************************************************************************************************

using blazemark::Benchmarks;
using blazemark::DynamicSparseRun;
using blazemark::Parser;




//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Type of a benchmark run.
//
// This type definition specifies the type of a single benchmark run for the benchmark for
// custom expression.
*/
using Run = DynamicSparseRun;
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
   const size_t N( run.getSize() );
   const size_t F( run.getNonZeros() );

   blaze::timing::WcTimer timer;
   double wct( 0.0 );
   size_t steps( 1UL );

   while( true ) {
      timer.start();
      blazemark::blaze::custom( N, F, steps );
      timer.end();
      wct = timer.last();
      if( wct >= 0.2 ) break;
      steps *= 2UL;
   }

   const size_t estimatedSteps( ( blazemark::runtime * steps ) / timer.last() );
   run.setSteps( blaze::max( 1UL, estimatedSteps ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  BENCHMARK FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Custom expression benchmark function.
//
// \param runs The specified benchmark runs.
// \param benchmarks The selection of benchmarks.
// \return void
*/
void custom( std::vector<Run>& runs, Benchmarks benchmarks )
{
   std::cout << std::left;

   std::sort( runs.begin(), runs.end() );

   size_t slowSize( blaze::inf );
   for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
      if( run->getSteps() == 0UL ) {
         if( run->getSize() < slowSize ) {
            estimateSteps( *run );
            if( run->getSteps() == 1UL )
               slowSize = run->getSize();
         }
         else run->setSteps( 1UL );
      }
   }

   if( benchmarks.runBlaze ) {
      std::cout << "   Blaze (Seconds):\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()     );
         const size_t F    ( run->getNonZeros() );
         const size_t steps( run->getSteps()    );
         run->setBlazeResult( blazemark::blaze::custom( N, F, steps ) );
         std::cout << "     " << std::setw(12) << run->getSize() << run->getBlazeResult() << std::endl;
      }
   }

#if BLAZEMARK_BOOST_MODE
   if( benchmarks.runBoost ) {
      std::cout << "   Boost uBLAS (Seconds):\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()     );
         const size_t F    ( run->getNonZeros() );
         const size_t steps( run->getSteps()    );
         run->setBoostResult( blazemark::boost::custom( N, F, steps ) );
         std::cout << "     " << std::setw(12) << run->getSize() << run->getBoostResult() << std::endl;
      }
   }
#endif

#if BLAZEMARK_BLITZ_MODE
   if( benchmarks.runBlitz ) {
      std::cout << "   Blitz++ (Seconds):\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()     );
         const size_t F    ( run->getNonZeros() );
         const size_t steps( run->getSteps()    );
         run->setBlitzResult( blazemark::blitz::custom( N, F, steps ) );
         std::cout << "     " << std::setw(12) << run->getSize() << run->getBlitzResult() << std::endl;
      }
   }
#endif

#if BLAZEMARK_GMM_MODE
   if( benchmarks.runGMM ) {
      std::cout << "   GMM++ (Seconds):\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()     );
         const size_t F    ( run->getNonZeros() );
         const size_t steps( run->getSteps()    );
         run->setGMMResult( blazemark::gmm::custom( N, F, steps ) );
         std::cout << "     " << std::setw(12) << run->getSize() << run->getGMMResult() << std::endl;
      }
   }
#endif

#if BLAZEMARK_ARMADILLO_MODE
   if( benchmarks.runArmadillo ) {
      std::cout << "   Armadillo (Seconds):\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()     );
         const size_t F    ( run->getNonZeros() );
         const size_t steps( run->getSteps()    );
         run->setArmadilloResult( blazemark::armadillo::custom( N, F, steps ) );
         std::cout << "     " << std::setw(12) << run->getSize() << run->getArmadilloResult() << std::endl;
      }
   }
#endif

#if BLAZEMARK_MTL_MODE
   if( benchmarks.runMTL ) {
      std::cout << "   MTL (Seconds):\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()     );
         const size_t F    ( run->getNonZeros() );
         const size_t steps( run->getSteps()    );
         run->setMTLResult( blazemark::mtl::custom( N, F, steps ) );
         std::cout << "     " << std::setw(12) << run->getSize() << run->getMTLResult() << std::endl;
      }
   }
#endif

#if BLAZEMARK_EIGEN_MODE
   if( benchmarks.runEigen ) {
      std::cout << "   Eigen (Seconds):\n";
      for( std::vector<Run>::iterator run=runs.begin(); run!=runs.end(); ++run ) {
         const size_t N    ( run->getSize()     );
         const size_t F    ( run->getNonZeros() );
         const size_t steps( run->getSteps()    );
         run->setEigenResult( blazemark::eigen::custom( N, F, steps ) );
         std::cout << "     " << std::setw(12) << run->getSize() << run->getEigenResult() << std::endl;
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
/*!\brief The main function for the benchmark for the custom expression.
//
// \param argc The total number of command line arguments.
// \param argv The array of command line arguments.
// \return void
*/
int main( int argc, char** argv )
{
   std::cout << "\n Custom Expression:\n";

   Benchmarks benchmarks;

   try {
      parseCommandLineArguments( argc, argv, benchmarks );
   }
   catch( std::exception& ex ) {
      std::cerr << "   " << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   const std::string installPath( INSTALL_PATH );
   const std::string parameterFile( installPath + "/params/custom.prm" );
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
      custom( runs, benchmarks );
   }
   catch( std::exception& ex ) {
      std::cerr << "   Error during benchmark execution: " << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
