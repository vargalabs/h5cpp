//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/IdentityMatrix.h
//  \brief Specialization of the Creator class template for IdentityMatrix
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_IDENTITYMATRIX_H_
#define _BLAZETEST_MATHTEST_CREATOR_IDENTITYMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/IdentityMatrix.h>
#include <blazetest/system/Types.h>


namespace blazetest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the Creator class template for identity \f$ N \times N \f$ matrices.
//
// This specialization of the Creator class template is able to create random identity matrices.
*/
template< typename T  // Element type of the identity matrix
        , bool SO >   // Storage order of the identity matrix
class Creator< blaze::IdentityMatrix<T,SO> >
{
 public:
   //**Type definitions****************************************************************************
   using Type = blaze::IdentityMatrix<T,SO>;  //!< Type to be created by the Creator.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Creator();
   explicit inline Creator( size_t n );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   // No explicitly declared copy assignment operator.

   blaze::IdentityMatrix<T,SO> operator()() const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t n_;  //!< The number of rows and columns of the identity matrix.
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
/*!\brief Constructor for the creator specialization for IdentityMatrix.
*/
template< typename T  // Element type of the identity matrix
        , bool SO >   // Storage order of the identity matrix
inline Creator< blaze::IdentityMatrix<T,SO> >::Creator()
   : n_( 3UL )  // The number of rows and columns of the identity matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the creator specialization for IdentityMatrix.
//
// \param n The number of rows and columns of the identity matrix.
*/
template< typename T  // Element type of the identity matrix
        , bool SO >   // Storage order of the identity matrix
inline Creator< blaze::IdentityMatrix<T,SO> >::Creator( size_t n )
   : n_( n )  // The number of rows and columns of the identity matrix
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created identity matrix.
//
// \return The randomly generated identity matrix.
*/
template< typename T  // Element type of the identity matrix
        , bool SO >   // Storage order of the identity matrix
inline blaze::IdentityMatrix<T,SO>
   Creator< blaze::IdentityMatrix<T,SO> >::operator()() const
{
   return blaze::IdentityMatrix<T,SO>( n_ );
}
//*************************************************************************************************

} // namespace blazetest

#endif
