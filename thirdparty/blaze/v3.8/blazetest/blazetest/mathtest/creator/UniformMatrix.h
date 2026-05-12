//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/UniformMatrix.h
//  \brief Specialization of the Creator class template for UniformMatrix
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_UNIFORMMATRIX_H_
#define _BLAZETEST_MATHTEST_CREATOR_UNIFORMMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/UniformMatrix.h>
#include <blazetest/mathtest/creator/Default.h>
#include <blazetest/mathtest/creator/Policies.h>
#include <blazetest/system/Types.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the Creator class template for uniform matrices.
//
// This specialization of the Creator class template is able to create random uniform matrices.
*/
template< typename T  // Element type of the uniform matrix
        , bool SO >   // Storage order of the uniform matrix
class Creator< blaze::UniformMatrix<T,SO> >
{
 public:
   //**Type definitions****************************************************************************
   using Type = blaze::UniformMatrix<T,SO>;  //!< Type to be created by the Creator.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Creator( const Creator<T>& elementCreator = Creator<T>() );
   explicit inline Creator( size_t m, size_t n, const Creator<T>& elementCreator = Creator<T>() );
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   // No explicitly declared copy assignment operator.

   blaze::UniformMatrix<T,SO> operator()() const;

   template< typename CP >
   blaze::UniformMatrix<T,SO> operator()( const CP& policy ) const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;       //!< The number of rows of the uniform matrix.
   size_t n_;       //!< The number of columns of the uniform matrix.
   Creator<T> ec_;  //!< Creator for the elements of the uniform matrix.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the creator specialization for UniformMatrix.
//
// \param elementCreator The creator for the elements of the uniform matrix.
*/
template< typename T  // Element type of the uniform matrix
        , bool SO >   // Storage order of the uniform matrix
inline Creator< blaze::UniformMatrix<T,SO> >::Creator( const Creator<T>& elementCreator )
   : m_( 3UL )              // The number of rows of the uniform matrix
   , n_( 3UL )              // The number of columns of the uniform matrix
   , ec_( elementCreator )  // Creator for the elements of the uniform matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the creator specialization for UniformMatrix.
//
// \param m The number of rows of the uniform matrix.
// \param n The number of columns of the uniform matrix.
// \param elementCreator The creator for the elements of the uniform matrix.
*/
template< typename T  // Element type of the uniform matrix
        , bool SO >   // Storage order of the uniform matrix
inline Creator< blaze::UniformMatrix<T,SO> >::Creator( size_t m, size_t n, const Creator<T>& elementCreator )
   : m_( m )                // The number of rows of the uniform matrix
   , n_( n )                // The number of columns of the uniform matrix
   , ec_( elementCreator )  // Creator for the elements of the uniform matrix
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created uniform matrix.
//
// \return The randomly generated uniform matrix.
*/
template< typename T  // Element type of the uniform matrix
        , bool SO >   // Storage order of the uniform matrix
inline blaze::UniformMatrix<T,SO> Creator< blaze::UniformMatrix<T,SO> >::operator()() const
{
   return (*this)( Default() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a randomly created uniform matrix.
//
// \param policy The creation policy for the elements of fundamental data type.
// \return The randomly generated uniform matrix.
*/
template< typename T     // Element type of the uniform matrix
        , bool SO >      // Storage order of the uniform matrix
template< typename CP >  // Creation policy
inline blaze::UniformMatrix<T,SO>
   Creator< blaze::UniformMatrix<T,SO> >::operator()( const CP& policy ) const
{
   return blaze::UniformMatrix<T,SO>( m_, n_, ec_( policy ) );
}
//*************************************************************************************************

} // namespace blazetest

#endif
