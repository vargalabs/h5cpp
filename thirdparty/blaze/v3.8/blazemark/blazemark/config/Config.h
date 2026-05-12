//=================================================================================================
/*!
//  \file blazemark/config/Config.h
//  \brief General configuration file for the blaze benchmark suite
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


//=================================================================================================
//
//  BENCHMARK CONFIGURATION
//
//=================================================================================================

//*************************************************************************************************
const bool runClike     ( true );  //!< C-like, manual implementation
const bool runClassic   ( true );  //!< Classical operator overloading benchmark tests
const bool runBLAS      ( true );  //!< BLAS benchmark tests
const bool runBlaze     ( true );  //!< Blaze benchmark tests
const bool runBoost     ( true );  //!< Boost uBLAS benchmark tests
const bool runBlitz     ( true );  //!< Blitz++ benchmark tests
const bool runGMM       ( true );  //!< GMM++ benchmark tests
const bool runArmadillo ( true );  //!< Armadillo benchmark tests
const bool runFLENS     ( true );  //!< FLENS benchmark tests
const bool runMTL       ( true );  //!< MTL benchmark tests
const bool runEigen     ( true );  //!< Eigen benchmark tests
//*************************************************************************************************


//*************************************************************************************************
const size_t reps     ( 3     );  //!< Configuration of the number of benchmark repetitions
const double runtime  ( 2.0   );  //!< Target runtime for a benchmark measurement
const double maxtime  ( 600.0 );  //!< Maximum runtime of a single benchmark measurement [s]
//*************************************************************************************************


//*************************************************************************************************
const double deviation( 5.0 );  //!< Maximum allowed deviation of the average benchmark time from the minimum time
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Element type for all benchmarks.
//
// This type definition specifies the data type used as element type for all dense and sparse
// vectors and matrices in all benchmarks. Possible settings for this type are <a>char, signed
// char, unsigned char, wchar_t, short, unsigned short, int, unsigned int, long, unsigned long,
// float, double, long double and complex</a> with the previously named types as element type.
// The default setting is \a double.
*/
using element_t = double;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Structure of sparse matrices
//
// This flag value specifies the structure of sparse matrices in all benchmarks. The structure
// can be either specified as blazemark::band, which results in the setup of banded matrices,
// or blazemark::random, which results in the setup of sparse matrices with randomly placed
// non-zero entries.
//
// Valid settings for the structure are blazemark::band and blaze::random.
*/
const MatrixStructure structure( random );
//*************************************************************************************************





//=================================================================================================
//
//  RANDOM NUMBER CONFIGURATION
//
//=================================================================================================

//*************************************************************************************************
const size_t seed( 128753984 );  //!< Seed for the random number generator.
//*************************************************************************************************




//=================================================================================================
//
//  EIGEN CONFIGURATION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Index data type for the Eigen sparse matrices.
//
// This type definition offers the possibility to specify the index type for the sparse matrices
// in all Eigen benchmarks.
*/
using EigenSparseIndexType = int;
//*************************************************************************************************
