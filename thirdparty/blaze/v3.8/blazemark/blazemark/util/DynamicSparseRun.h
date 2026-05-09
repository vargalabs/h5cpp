//=================================================================================================
/*!
//  \file blazemark/util/DynamicSparseRun.h
//  \brief Header file for the DynamicSparseRun class
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

#ifndef _BLAZEMARK_UTIL_DYNAMICSPARSERUN_H_
#define _BLAZEMARK_UTIL_DYNAMICSPARSERUN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iomanip>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <blaze/math/Infinity.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/UnsignedValue.h>
#include <blazemark/system/Types.h>


namespace blazemark {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Data structure for the parameters of a benchmark run with sparse vectors and/or matrices.
//
// This auxiliary data structure represents the necessary parameters for a benchmark run with
// sparse vectors and/or matrices.
*/
class DynamicSparseRun
{
 private:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline DynamicSparseRun();
   //@}
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline DynamicSparseRun( size_t size, size_t nonzeros );
   explicit inline DynamicSparseRun( size_t size, size_t nonzeros, size_t steps );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Copy assignment operator********************************************************************
   // No explicitly declared copy assignment operator.
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t getSize           () const;
   inline size_t getNonZeros       () const;
   inline float  getFillingDegree  () const;
   inline size_t getSteps          () const;
   inline size_t getFlops          () const;
   inline double getClikeResult    () const;
   inline double getClassicResult  () const;
   inline double getBlazeResult    () const;
   inline double getBoostResult    () const;
   inline double getBlitzResult    () const;
   inline double getGMMResult      () const;
   inline double getArmadilloResult() const;
   inline double getFLENSResult    () const;
   inline double getMTLResult      () const;
   inline double getEigenResult    () const;

   inline void   setSize    ( size_t newSize     );
   inline void   setNonZeros( size_t newNonZeros );
   inline void   setSteps   ( size_t newSteps    );
   inline void   setFlops   ( size_t newFlops    );
   inline void   setClikeResult    ( double result );
   inline void   setClassicResult  ( double result );
   inline void   setBlazeResult    ( double result );
   inline void   setBoostResult    ( double result );
   inline void   setBlitzResult    ( double result );
   inline void   setGMMResult      ( double result );
   inline void   setArmadilloResult( double result );
   inline void   setFLENSResult    ( double result );
   inline void   setMTLResult      ( double result );
   inline void   setEigenResult    ( double result );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t size_;       //!< The target size of the sparse vectors/matrices.
                       /*!< In case of a sparse vector, this value directly corresponds to the size
                            of the vector, in case of a sparse matrix \a size_ corresponds to the
                            number of rows and columns. */
   size_t nonzeros_;   //!< The number of non-zero elements in the sparse vectors/matrices.
                       /*!< In case of a sparse vector, this value directly corresponds to the
                            number of sparse elements. In case of a sparse matrix, the value
                            specifies the number of non-zero elements per row. */
   size_t steps_;      //!< The number of steps for the benchmark run.
                       /*!< The (composite) arithmetic operation of each benchmark is run several
                            times to guarantee reasonable runtimes. \a steps_ corresponds to the
                            number of performed iterations. */
   size_t flops_;      //!< The number of flops required for the benchmark run.
                       /*!< This value corresponds to the total number of floating point operations
                            (Flops) required for a single computation of the (composite) arithmetic
                            operation. */
   double clike_;      //!< Benchmark result of the C-like implementation.
   double classic_;    //!< Benchmark result of classic C++ operator overloading.
   double blaze_;      //!< Benchmark result of the Blaze library.
   double boost_;      //!< Benchmark result of the Boost uBLAS library.
   double blitz_;      //!< Benchmark result of the Blitz++ library.
   double gmm_;        //!< Benchmark result of the GMM++ library.
   double armadillo_;  //!< Benchmark result of the Armadillo library.
   double flens_;      //!< Benchmark result of the FLENS library.
   double mtl_;        //!< Benchmark result of the MTL4 library.
   double eigen_;      //!< Benchmark result of the Eigen3 library.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename > friend class Parser;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default constructor for the DynamicSparseRun class.
//
// The default constructor in exclusively accessible for the blazemark::Parser class.
*/
inline DynamicSparseRun::DynamicSparseRun()
   : size_     ( 0UL )  // The target size of the sparse vectors/matrices
   , nonzeros_ ( 0UL )  // The number of non-zero elements in the sparse vectors/matrices.
   , steps_    ( 0UL )  // The number of steps for the benchmark run
   , flops_    ( 0UL )  // The number of flops required for the benchmark run
   , clike_    ( 0.0 )  // Benchmark result of the C-like implementation
   , classic_  ( 0.0 )  // Benchmark result of the classic C++ implementation
   , blaze_    ( 0.0 )  // Benchmark result of the Blaze library
   , boost_    ( 0.0 )  // Benchmark result of the Boost uBLAS library
   , blitz_    ( 0.0 )  // Benchmark result of the Blitz++ library
   , gmm_      ( 0.0 )  // Benchmark result of the GMM++ library
   , armadillo_( 0.0 )  // Benchmark result of the Armadillo library
   , flens_    ( 0.0 )  // Benchmark result of the FLENS library
   , mtl_      ( 0.0 )  // Benchmark result of the MTL4 library
   , eigen_    ( 0.0 )  // Benchmark result of the Eigen3 library
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Two-argument constructor for the DynamicSparseRun class.
//
// \param size The size of the sparse vector or matrix \f$ [1..\infty) \f$.
// \param nonzeros The number of non-zero elements in the sparse vectors/matrices \f$ [1..size] \f$.
// \exception std::invalid_argument Invalid size parameter.
// \exception std::invalid_argument Invalid number of non-zero elements.
//
// This constructor creates a sparse run with a specified target size for the sparse vectors
// and/or matrices and a specified number of non-zero elements. The number of steps will
// automatically be evaluated to guarantuee an approximate runtime of blazemark::runtime
// seconds (see for the 'blazemark/config/Config.h' file for more details).
*/
inline DynamicSparseRun::DynamicSparseRun( size_t size, size_t nonzeros )
   : size_     ( size     )  // The target size of the sparse vectors/matrices
   , nonzeros_ ( nonzeros )  // The number of non-zero elements in the sparse vectors/matrices.
   , steps_    ( 0UL      )  // The number of steps for the benchmark run
   , flops_    ( 0UL      )  // The number of flops required for the benchmark run
   , clike_    ( 0.0      )  // Benchmark result of the C-like implementation
   , classic_  ( 0.0      )  // Benchmark result of the classic C++ implementation
   , blaze_    ( 0.0      )  // Benchmark result of the Blaze library
   , boost_    ( 0.0      )  // Benchmark result of the Boost uBLAS library
   , blitz_    ( 0.0      )  // Benchmark result of the Blitz++ library
   , gmm_      ( 0.0      )  // Benchmark result of the GMM++ library
   , armadillo_( 0.0      )  // Benchmark result of the Armadillo library
   , flens_    ( 0.0      )  // Benchmark result of the FLENS library
   , mtl_      ( 0.0      )  // Benchmark result of the MTL4 library
   , eigen_    ( 0.0      )  // Benchmark result of the Eigen3 library
{
   // Checking the target size for the sparse vectors/matrices
   if( size_ == size_t(0) )
      throw std::invalid_argument( "Invalid size parameter" );

   // Checking the number of non-zero elements
   if( nonzeros_ == size_t(0) || nonzeros_ > size_ )
      throw std::invalid_argument( "Invalid number of non-zero elements" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Three-argument constructor for the DynamicSparseRun class.
//
// \param size The size of the sparse vector or matrix \f$ [1..\infty) \f$.
// \param nonzeros The number of non-zero elements in the sparse vectors/matrices \f$ [1..size] \f$.
// \param steps The number of steps for the benchmark \f$ [1..\infty) \f$.
// \exception std::invalid_argument Invalid size parameter.
// \exception std::invalid_argument Invalid number of non-zero elements.
//
// This constructor creates a sparse run with a specified target size for the sparse vectors
// and/or matrcies, a specfied number of non-zero elements, and a specified number of steps
// for the benchmark. In case \a steps is zero, the number of steps will automatically be
// evaluated to guarantee an approximate runtime of blazemark::runtime seconds (see for
// the 'blazemark/config/Config.h' file for more details).
*/
inline DynamicSparseRun::DynamicSparseRun( size_t size, size_t nonzeros, size_t steps )
   : size_     ( size     )  // The target size of the sparse vectors/matrices
   , nonzeros_ ( nonzeros )  // The number of non-zero elements in the sparse vectors/matrices
   , steps_    ( steps    )  // The number of steps for the benchmark run
   , flops_    ( 0UL      )  // The number of flops required for the benchmark run
   , clike_    ( 0.0      )  // Benchmark result of the C-like implementation
   , classic_  ( 0.0      )  // Benchmark result of the classic C++ implementation
   , blaze_    ( 0.0      )  // Benchmark result of the Blaze library
   , boost_    ( 0.0      )  // Benchmark result of the Boost uBLAS library
   , blitz_    ( 0.0      )  // Benchmark result of the Blitz++ library
   , gmm_      ( 0.0      )  // Benchmark result of the GMM++ library
   , armadillo_( 0.0      )  // Benchmark result of the Armadillo library
   , flens_    ( 0.0      )  // Benchmark result of the FLENS library
   , mtl_      ( 0.0      )  // Benchmark result of the MTL4 library
   , eigen_    ( 0.0      )  // Benchmark result of the Eigen3 library
{
   // Checking the target size for the sparse vectors/matrices
   if( size_ == size_t(0) )
      throw std::invalid_argument( "Invalid size parameter" );

   // Checking the number of non-zero elements
   if( nonzeros_ == size_t(0) || nonzeros_ > size_ )
      throw std::invalid_argument( "Invalid number of non-zero elements" );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the target size of the sparse vectors/matrices of the benchmark run.
//
// \return The target size of the vectors/matrices.
*/
inline size_t DynamicSparseRun::getSize() const
{
   return size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements of the sparse vectors/matrices of the benchmark run.
//
// \return The number of non-zero elements of the vectors/matrices.
*/
inline size_t DynamicSparseRun::getNonZeros() const
{
   return nonzeros_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the filling degree of the sparse vectors/matrices of the benchmark run.
//
// \return The filling degree of the vectors/matrices in percent.
*/
inline float DynamicSparseRun::getFillingDegree() const
{
   return float(nonzeros_) / float(size_) * 100.0F;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of steps of the benchmark run.
//
// \return The number of steps of the benchmark run.
*/
inline size_t DynamicSparseRun::getSteps() const
{
   return steps_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of required floating point operations.
//
// \return The number of required floating point operations.
*/
inline size_t DynamicSparseRun::getFlops() const
{
   return flops_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the C-like implementation.
//
// \return The result of the C-like implementation.
*/
inline double DynamicSparseRun::getClikeResult() const
{
   return clike_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the classic C++ implementation.
//
// \return The result of the classic C++ implementation.
*/
inline double DynamicSparseRun::getClassicResult() const
{
   return classic_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the Blaze library.
//
// \return The result of the Blaze library.
*/
inline double DynamicSparseRun::getBlazeResult() const
{
   return blaze_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the Boost uBLAS library.
//
// \return The result of the Boost uBLAS library.
*/
inline double DynamicSparseRun::getBoostResult() const
{
   return boost_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the Blitz++ library.
//
// \return The result of the Blitz++ library.
*/
inline double DynamicSparseRun::getBlitzResult() const
{
   return blitz_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the GMM++ library.
//
// \return The result of the GMM++ library.
*/
inline double DynamicSparseRun::getGMMResult() const
{
   return gmm_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the Armadillo library.
//
// \return The result of the Armadillo library.
*/
inline double DynamicSparseRun::getArmadilloResult() const
{
   return armadillo_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the FLENS library.
//
// \return The result of the FLENS library.
*/
inline double DynamicSparseRun::getFLENSResult() const
{
   return flens_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the MTL4 library.
//
// \return The result of the MTL4 library.
*/
inline double DynamicSparseRun::getMTLResult() const
{
   return mtl_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the Eigen3 library.
//
// \return The result of the Eigen3 library.
*/
inline double DynamicSparseRun::getEigenResult() const
{
   return eigen_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the target size of the vectors/matrices of the benchmark run.
//
// \param newSize The new target size of the vectors/matrices of the benchmark run.
// \return void
// \exception std::invalid_argument Invalid size parameter.
*/
inline void DynamicSparseRun::setSize( size_t newSize )
{
   if( newSize == size_t(0) )
      throw std::invalid_argument( "Invalid size parameter" );
   size_ = newSize;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the number of non-zero elements of the vectors/matrices of the benchmark run.
//
// \param newNonZeros The new number of non-zero elements of the vectors/matrices of the benchmark run.
// \return void
// \exception std::invalid_argument Invalid number of non-zero elements.
*/
inline void DynamicSparseRun::setNonZeros( size_t newNonZeros )
{
   if( newNonZeros == size_t(0) || newNonZeros > size_ )
      throw std::invalid_argument( "Invalid number of non-zero elements" );
   nonzeros_ = newNonZeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the number of steps of the benchmark run.
//
// \param newSteps The new number of steps of the benchmark run.
// \return void
*/
inline void DynamicSparseRun::setSteps( size_t newSteps )
{
   steps_ = newSteps;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the number of required floating point operations.
//
// \param newFlops The new number of required floating point operations.
// \return void
*/
inline void DynamicSparseRun::setFlops( size_t newFlops )
{
   flops_ = newFlops;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the benchmark result of the C-like implementation.
//
// \param result The result of the C-like implementation.
// \return void
// \exception std::invalid_argument Invalid result value.
*/
inline void DynamicSparseRun::setClikeResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   clike_ = result;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the benchmark result of the classic C++ implementation.
//
// \param result The result of the classic C++ implementation.
// \return void
// \exception std::invalid_argument Invalid result value.
*/
inline void DynamicSparseRun::setClassicResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   classic_ = result;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the benchmark result of the Blaze library.
//
// \param result The result of the Blaze library.
// \return void
// \exception std::invalid_argument Invalid result value.
*/
inline void DynamicSparseRun::setBlazeResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   blaze_ = result;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the benchmark result of the Boost uBLAS library.
//
// \param result The result of the Boost uBLAS library.
// \return void
// \exception std::invalid_argument Invalid result value.
*/
inline void DynamicSparseRun::setBoostResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   boost_ = result;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the benchmark result of the Blitz++ library.
//
// \param result The result of the Blitz++ library.
// \return void
// \exception std::invalid_argument Invalid result value.
*/
inline void DynamicSparseRun::setBlitzResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   blitz_ = result;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the benchmark result of the GMM++ library.
//
// \param result The result of the GMM++ library.
// \return void
// \exception std::invalid_argument Invalid result value.
*/
inline void DynamicSparseRun::setGMMResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   gmm_ = result;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the benchmark result of the Armadillo library.
//
// \param result The result of the Armadillo library.
// \return void
// \exception std::invalid_argument Invalid result value.
*/
inline void DynamicSparseRun::setArmadilloResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   armadillo_ = result;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the benchmark result of the FLENS library.
//
// \param result The result of the FLENS library.
// \return void
// \exception std::invalid_argument Invalid result value.
*/
inline void DynamicSparseRun::setFLENSResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   flens_ = result;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the benchmark result of the MTL4 library.
//
// \param result The result of the MTL4 library.
// \return void
// \exception std::invalid_argument Invalid result value.
*/
inline void DynamicSparseRun::setMTLResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   mtl_ = result;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the benchmark result of the Eigen3 library.
//
// \param result The result of the Eigen3 library.
// \return void
// \exception std::invalid_argument Invalid result value.
*/
inline void DynamicSparseRun::setEigenResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   eigen_ = result;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Less-than comparison between two DynamicSparseRun objects.
//
// \param lhs The left-hand side DynamicSparseRun object.
// \param rhs The right-hand side DynamicSparseRun object.
// \return \a true if the left value is less than the right value, \a false if not.
//
// DynamicSparseRun objects are primarily sorted according to the number of non-zero elements
// and secondarily according to the size of the vectors/matrices: In case the nonzeros value
// of the left-hand side DynamicSparseRun object is smaller or in case both nonzeros values
// are equal and the size value of the left-hand side is smaller the function returns \a true.
// Otherwise the function returns \a false.
*/
inline bool operator<( const DynamicSparseRun& lhs, const DynamicSparseRun& rhs )
{
   return ( lhs.getFillingDegree() < rhs.getFillingDegree() ) ||
          ( ::blaze::equal( lhs.getFillingDegree(), rhs.getFillingDegree() ) && lhs.getSize() < rhs.getSize() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for the DynamicSparseRun class.
//
// \param os Reference to the output stream.
// \param run Reference to a DynamicSparseRun object.
// \return The output stream.
*/
inline std::ostream& operator<<( std::ostream& os, const DynamicSparseRun& run )
{
   const std::ios::fmtflags flags( os.flags() );

   os << std::left << "   N=" << run.getSize() << ", fill level=" << run.getFillingDegree() << "%, steps=" << run.getSteps() << "\n";

   const double clike    ( run.getClikeResult()     );
   const double classic  ( run.getClassicResult()   );
   const double blaze    ( run.getBlazeResult()     );
   const double boost    ( run.getBoostResult()     );
   const double blitz    ( run.getBlitzResult()     );
   const double gmm      ( run.getGMMResult()       );
   const double armadillo( run.getArmadilloResult() );
   const double flens    ( run.getFLENSResult()     );
   const double mtl      ( run.getMTLResult()       );
   const double eigen    ( run.getEigenResult()     );

   double minTime = ::blaze::inf;

   if( clike     != 0.0 ) minTime = ::blaze::min( minTime, clike     );
   if( classic   != 0.0 ) minTime = ::blaze::min( minTime, classic   );
   if( blaze     != 0.0 ) minTime = ::blaze::min( minTime, blaze     );
   if( boost     != 0.0 ) minTime = ::blaze::min( minTime, boost     );
   if( blitz     != 0.0 ) minTime = ::blaze::min( minTime, blitz     );
   if( gmm       != 0.0 ) minTime = ::blaze::min( minTime, gmm       );
   if( armadillo != 0.0 ) minTime = ::blaze::min( minTime, armadillo );
   if( flens     != 0.0 ) minTime = ::blaze::min( minTime, flens     );
   if( mtl       != 0.0 ) minTime = ::blaze::min( minTime, mtl       );
   if( eigen     != 0.0 ) minTime = ::blaze::min( minTime, eigen     );

   if( clike     != 0.0 ) os << "     C-like      = " << std::setw(8) << ( clike     / minTime ) << " (" << clike     << ")\n";
   if( classic   != 0.0 ) os << "     Classic     = " << std::setw(8) << ( classic   / minTime ) << " (" << classic   << ")\n";
   if( blaze     != 0.0 ) os << "     Blaze       = " << std::setw(8) << ( blaze     / minTime ) << " (" << blaze     << ")\n";
   if( boost     != 0.0 ) os << "     Boost uBLAS = " << std::setw(8) << ( boost     / minTime ) << " (" << boost     << ")\n";
   if( blitz     != 0.0 ) os << "     Blitz++     = " << std::setw(8) << ( blitz     / minTime ) << " (" << blitz     << ")\n";
   if( gmm       != 0.0 ) os << "     GMM++       = " << std::setw(8) << ( gmm       / minTime ) << " (" << gmm       << ")\n";
   if( armadillo != 0.0 ) os << "     Armadillo   = " << std::setw(8) << ( armadillo / minTime ) << " (" << armadillo << ")\n";
   if( flens     != 0.0 ) os << "     FLENS       = " << std::setw(8) << ( flens     / minTime ) << " (" << flens     << ")\n";
   if( mtl       != 0.0 ) os << "     MTL         = " << std::setw(8) << ( mtl       / minTime ) << " (" << mtl       << ")\n";
   if( eigen     != 0.0 ) os << "     Eigen       = " << std::setw(8) << ( eigen     / minTime ) << " (" << eigen     << ")\n";

   os << std::flush;

   os.flags( flags );
   return os;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global input operator for the DynamicSparseRun class.
//
// \param is Reference to the input stream.
// \param run Reference to a DynamicSparseRun object.
// \return The input stream.
//
// The input operator guarantees that this object is not changed in the case of an input error.
// Only values suitable for the according built-in unsigned integral data type \a T are allowed.
// Otherwise, the input stream's position is returned to its previous position and the
// \a std::istream::failbit is set.
*/
inline std::istream& operator>>( std::istream& is, DynamicSparseRun& run )
{
   char c1, c2, c3;
   ::blaze::UnsignedValue<size_t> size, nonzeros, steps;
   const std::istream::pos_type pos( is.tellg() );

   if( !(is >> c1 >> size >> c2 >> nonzeros >> c3) ||
       c1 != '(' || size == 0 || c2 != ',' || nonzeros == 0 ||
       ( c3 != ')' && ( c3 != ',' || ( !(is >> steps >> c1) || c1 != ')' || steps == 0 ) ) ) )
   {
      is.clear();
      is.seekg( pos );
      is.setstate( std::istream::failbit );
      return is;
   }

   run.setSize( size );
   run.setNonZeros( nonzeros );
   run.setSteps( steps );

   return is;
}
//*************************************************************************************************

} // namespace blazemark

#endif
