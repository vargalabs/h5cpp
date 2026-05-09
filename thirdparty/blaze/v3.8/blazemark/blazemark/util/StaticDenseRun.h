//=================================================================================================
/*!
//  \file blazemark/util/StaticDenseRun.h
//  \brief Header file for the StaticDenseRun class
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

#ifndef _BLAZEMARK_UTIL_STATICDENSERUN_H_
#define _BLAZEMARK_UTIL_STATICDENSERUN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iomanip>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <blaze/math/Infinity.h>
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
/*!\brief Data structure for the parameters of a benchmark run with fixed size dense vectors/matrices.
//
// This auxiliary data structure represents the necessary parameters for a benchmark run with
// fixed size dense vectors and/or matrices.
*/
template< size_t N >  // Fixed size of the vectors/matrices
class StaticDenseRun
{
 private:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline StaticDenseRun();
   //@}
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline StaticDenseRun( size_t number );
   explicit inline StaticDenseRun( size_t number, size_t steps );
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
   inline size_t getSize  () const;
   inline size_t getNumber() const;
   inline size_t getSteps () const;
   inline size_t getFlops () const;
   inline double getClikeResult    () const;
   inline double getClassicResult  () const;
   inline double getBLASResult     () const;
   inline double getBlazeResult    () const;
   inline double getBoostResult    () const;
   inline double getBlitzResult    () const;
   inline double getGMMResult      () const;
   inline double getArmadilloResult() const;
   inline double getFLENSResult    () const;
   inline double getMTLResult      () const;
   inline double getEigenResult    () const;

   inline void   setNumber( size_t newNumber );
   inline void   setSteps ( size_t newSteps  );
   inline void   setFlops ( size_t newFlops  );
   inline void   setClikeResult    ( double result );
   inline void   setClassicResult  ( double result );
   inline void   setBLASResult     ( double result );
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
   size_t number_;     //!< The target number of fixed size vectors/matrices.
                       /*!< This value specifies the number of (potentially very small sized)
                            vectors and/or matrices used in the benchmark. */
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
   double blas_;       //!< Benchmark result of the BLAS implementation.
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
/*!\brief Default constructor for the StaticDenseRun class.
//
// The default constructor in exclusively accessible for the blazemark::Parser class.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline StaticDenseRun<N>::StaticDenseRun()
   : number_   ( 0UL )  // The target number of fixed size vectors/matrices
   , steps_    ( 0UL )  // The number of steps for the benchmark run
   , flops_    ( 0UL )  // The number of flops required for the benchmark run
   , clike_    ( 0.0 )  // Benchmark result of the C-like implementation
   , classic_  ( 0.0 )  // Benchmark result of the classic C++ implementation
   , blas_     ( 0.0 )  // Benchmark result of the BLAS implementation
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
/*!\brief Single-argument constructor for the StaticDenseRun class.
//
// \param number The number of fixed size vectors and/or matrices \f$ [1..\infty) \f$.
// \exception std::invalid_argument Invalid number parameter.
//
// This constructor creates a dense run with a specified number of fixed size vectors and/or
// matrices. The number of steps will automatically be evaluated to guarantuee an approximate
// runtime of blazemark::runtime seconds (see for the 'blazemark/config/Config.h' file for
// more details).
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline StaticDenseRun<N>::StaticDenseRun( size_t number )
   : number_   ( number )  // The target number of fixed size vectors/matrices
   , steps_    ( 0UL    )  // The number of steps for the benchmark run
   , flops_    ( 0UL    )  // The number of flops required for the benchmark run
   , clike_    ( 0.0    )  // Benchmark result of the C-like implementation
   , classic_  ( 0.0    )  // Benchmark result of the classic C++ implementation
   , blas_     ( 0.0    )  // Benchmark result of the BLAS implementation
   , blaze_    ( 0.0    )  // Benchmark result of the Blaze library
   , boost_    ( 0.0    )  // Benchmark result of the Boost uBLAS library
   , blitz_    ( 0.0    )  // Benchmark result of the Blitz++ library
   , gmm_      ( 0.0    )  // Benchmark result of the GMM++ library
   , armadillo_( 0.0    )  // Benchmark result of the Armadillo library
   , flens_    ( 0.0    )  // Benchmark result of the FLENS library
   , mtl_      ( 0.0    )  // Benchmark result of the MTL4 library
   , eigen_    ( 0.0    )  // Benchmark result of the Eigen3 library
{
   // Checking the target number for the fixed size vectors/matrices
   if( number_ == size_t(0) )
      throw std::invalid_argument( "Invalid number parameter" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Two-argument constructor for the StaticDenseRun class.
//
// \param number The number of fixed size vectors and/or matrices \f$ [1..\infty) \f$.
// \param steps The number of steps for the benchmark \f$ [1..\infty) \f$.
// \exception std::invalid_argument Invalid number parameter.
//
// This constructor creates a dense run with a specified number of fixed size vectors and/or
// matrices and a specified number of steps for the benchmark. In case \a steps is zero, the
// number of steps will automatically be evaluated to guarantee an approximate runtime of
// blazemark::runtime seconds (see for the 'blazemark/config/Config.h' file for more details).
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline StaticDenseRun<N>::StaticDenseRun( size_t number, size_t steps )
   : number_   ( number )  // The target number of fixed size vectors/matrices
   , steps_    ( steps  )  // The number of steps for the benchmark run
   , flops_    ( 0UL    )  // The number of flops required for the benchmark run
   , clike_    ( 0.0    )  // Benchmark result of the C-like implementation
   , classic_  ( 0.0    )  // Benchmark result of the classic C++ implementation
   , blas_     ( 0.0    )  // Benchmark result of the BLAS implementation
   , blaze_    ( 0.0    )  // Benchmark result of the Blaze library
   , boost_    ( 0.0    )  // Benchmark result of the Boost uBLAS library
   , blitz_    ( 0.0    )  // Benchmark result of the Blitz++ library
   , gmm_      ( 0.0    )  // Benchmark result of the GMM++ library
   , armadillo_( 0.0    )  // Benchmark result of the Armadillo library
   , flens_    ( 0.0    )  // Benchmark result of the FLENS library
   , mtl_      ( 0.0    )  // Benchmark result of the MTL4 library
   , eigen_    ( 0.0    )  // Benchmark result of the Eigen3 library
{
   // Checking the target number for the fixed size vectors/matrices
   if( number_ == size_t(0) )
      throw std::invalid_argument( "Invalid number parameter" );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the target size of the dense vectors/matrices of the benchmark run.
//
// \return The target size of the vectors/matrices.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline size_t StaticDenseRun<N>::getSize() const
{
   return N;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of fixed size vectors/matrices of the benchmark run.
//
// \return The number of fixed size vectors/matrices.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline size_t StaticDenseRun<N>::getNumber() const
{
   return number_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of steps of the benchmark run.
//
// \return The number of steps of the benchmark run.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline size_t StaticDenseRun<N>::getSteps() const
{
   return steps_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of required floating point operations.
//
// \return The number of required floating point operations.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline size_t StaticDenseRun<N>::getFlops() const
{
   return flops_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the C-like implementation.
//
// \return The result of the C-like implementation.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline double StaticDenseRun<N>::getClikeResult() const
{
   return clike_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the classic C++ implementation.
//
// \return The result of the classic C++ implementation.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline double StaticDenseRun<N>::getClassicResult() const
{
   return classic_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the BLAS implementation.
//
// \return The result of the BLAS implementation.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline double StaticDenseRun<N>::getBLASResult() const
{
   return blas_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the Blaze library.
//
// \return The result of the Blaze library.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline double StaticDenseRun<N>::getBlazeResult() const
{
   return blaze_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the Boost uBLAS library.
//
// \return The result of the Boost uBLAS library.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline double StaticDenseRun<N>::getBoostResult() const
{
   return boost_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the Blitz++ library.
//
// \return The result of the Blitz++ library.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline double StaticDenseRun<N>::getBlitzResult() const
{
   return blitz_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the GMM++ library.
//
// \return The result of the GMM++ library.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline double StaticDenseRun<N>::getGMMResult() const
{
   return gmm_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the Armadillo library.
//
// \return The result of the Armadillo library.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline double StaticDenseRun<N>::getArmadilloResult() const
{
   return armadillo_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the FLENS library.
//
// \return The result of the FLENS library.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline double StaticDenseRun<N>::getFLENSResult() const
{
   return flens_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the MTL4 library.
//
// \return The result of the MTL4 library.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline double StaticDenseRun<N>::getMTLResult() const
{
   return mtl_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the benchmark result of the Eigen3 library.
//
// \return The result of the Eigen3 library.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline double StaticDenseRun<N>::getEigenResult() const
{
   return eigen_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the number of fixed size vectors/matrices of the benchmark run.
//
// \param newNumber The new number of fixed size vectors/matrices of the benchmark run.
// \return void
// \exception std::invalid_argument Invalid number parameter.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setNumber( size_t newNumber )
{
   if( newNumber == size_t(0) )
      throw std::invalid_argument( "Invalid number parameter" );
   number_ = newNumber;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the number of steps of the benchmark run.
//
// \param newSteps The new number of steps of the benchmark run.
// \return void
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setSteps( size_t newSteps )
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
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setFlops( size_t newFlops )
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
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setClikeResult( double result )
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
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setClassicResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   classic_ = result;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the benchmark result of the BLAS implementation.
//
// \param result The result of the BLAS implementation.
// \return void
// \exception std::invalid_argument Invalid result value.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setBLASResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   blas_ = result;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the benchmark result of the Blaze library.
//
// \param result The result of the Blaze library.
// \return void
// \exception std::invalid_argument Invalid result value.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setBlazeResult( double result )
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
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setBoostResult( double result )
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
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setBlitzResult( double result )
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
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setGMMResult( double result )
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
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setArmadilloResult( double result )
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
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setFLENSResult( double result )
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
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setMTLResult( double result )
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
template< size_t N >  // Fixed size of the vectors/matrices
inline void StaticDenseRun<N>::setEigenResult( double result )
{
   if( result < 0.0 )
      throw std::invalid_argument( "Invalid result value" );
   eigen_ = result;
}
//*************************************************************************************************




//=================================================================================================
//
//  UNDEFINED CLASS TEMPLATE SPECIALIZATION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of StaticDenseRun for 0 elements.
//
// This specialization of the StaticDenseRun class template is left undefined and therefore
// prevents the instantiation for 0 elements.
*/
template<>
class StaticDenseRun<0UL>;
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Less-than comparison between two StaticDenseRun objects.
//
// \param lhs The left-hand side StaticDenseRun object.
// \param rhs The right-hand side StaticDenseRun object.
// \return \a true if the left value is less than the right value, \a false if not.
//
// StaticDenseRun objects are sorted according to the size value: In case the size value of
// the left-hand side StaticDenseRun object is smaller than the size value of the right-hand
// size StaticDenseRun object, the function returns \a true. Otherwise it returns \a false.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline bool operator<( const StaticDenseRun<N>& lhs, const StaticDenseRun<N>& rhs )
{
   return lhs.getSize() < rhs.getSize();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for the StaticDenseRun class.
//
// \param os Reference to the output stream.
// \param run Reference to a StaticDenseRun object.
// \return The output stream.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline std::ostream& operator<<( std::ostream& os, const StaticDenseRun<N>& run )
{
   const std::ios::fmtflags flags( os.flags() );

   os << std::left << "   N=" << run.getSize() << ", number=" << run.getNumber() << ", steps=" << run.getSteps() << "\n";

   const double clike    ( run.getClikeResult()     );
   const double classic  ( run.getClassicResult()   );
   const double blas     ( run.getBLASResult()      );
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
   if( blas      != 0.0 ) minTime = ::blaze::min( minTime, blas      );
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
   if( blas      != 0.0 ) os << "     BLAS        = " << std::setw(8) << ( blas      / minTime ) << " (" << blas      << ")\n";
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
/*!\brief Global input operator for the StaticDenseRun class.
//
// \param is Reference to the input stream.
// \param run Reference to a StaticDenseRun object.
// \return The input stream.
//
// The input operator guarantees that this object is not changed in the case of an input error.
// Only values suitable for the according built-in unsigned integral data type \a T are allowed.
// Otherwise, the input stream's position is returned to its previous position and the
// \a std::istream::failbit is set.
*/
template< size_t N >  // Fixed size of the vectors/matrices
inline std::istream& operator>>( std::istream& is, StaticDenseRun<N>& run )
{
   char c1, c2;
   ::blaze::UnsignedValue<size_t> number, steps;
   const std::istream::pos_type pos( is.tellg() );

   if( !(is >> c1 >> number >> c2) || c1 != '(' || number == 0 ||
       ( c2 != ')' && ( c2 != ',' || ( !(is >> steps >> c1) || c1 != ')' || steps == 0 ) ) ) )
   {
      is.clear();
      is.seekg( pos );
      is.setstate( std::istream::failbit );
      return is;
   }

   run.setNumber( number );
   run.setSteps ( steps  );

   return is;
}
//*************************************************************************************************

} // namespace blazemark

#endif
