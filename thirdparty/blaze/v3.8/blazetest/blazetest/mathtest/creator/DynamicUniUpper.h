//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/DynamicUniUpper.h
//  \brief Specialization of the Creator class template for UniUpperMatrix<DynamicMatrix>
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_DYNAMICUNIUPPER_H_
#define _BLAZETEST_MATHTEST_CREATOR_DYNAMICUNIUPPER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
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
/*!\brief Specialization of the Creator class template for uniupper dynamic matrices.
//
// This specialization of the Creator class template is able to create random uniupper dynamic
// matrices.
*/
template< typename T  // Element type of the dynamic matrix
        , bool SO >   // Storage order of the dynamic matrix
class Creator< blaze::UniUpperMatrix< blaze::DynamicMatrix<T,SO> > >
{
 public:
   //**Type definitions****************************************************************************
   //! Type to be created by the Creator.
   using Type = blaze::UniUpperMatrix< blaze::DynamicMatrix<T,SO> >;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Creator( const Creator<T>& elementCreator = Creator<T>() );
   explicit inline Creator( size_t n, const Creator<T>& elementCreator = Creator<T>() );
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

   blaze::UniUpperMatrix< blaze::DynamicMatrix<T,SO> > operator()() const;

   template< typename CP >
   blaze::UniUpperMatrix< blaze::DynamicMatrix<T,SO> > operator()( const CP& policy ) const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t n_;       //!< The number of rows and columns of the uniupper dynamic matrix.
   Creator<T> ec_;  //!< Creator for the elements of the uniupper dynamic matrix.
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
/*!\brief Constructor for the creator specialization for UniUpperMatrix<DynamicMatrix>.
//
// \param elementCreator The creator for the elements of the uniupper dynamic matrix.
*/
template< typename T  // Element type of the dynamic matrix
        , bool SO >   // Storage order of the dynamic matrix
inline Creator< blaze::UniUpperMatrix< blaze::DynamicMatrix<T,SO> > >::Creator( const Creator<T>& elementCreator )
   : n_( 3UL )              // The number of rows and columns of the uniupper dynamic matrix
   , ec_( elementCreator )  // Creator for the elements of the uniupper dynamic matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the creator specialization for UniUpperMatrix<DynamicMatrix>.
//
// \param n The number of rows and columns of the uniupper dynamic matrix.
// \param elementCreator The creator for the elements of the uniupper dynamic matrix.
*/
template< typename T  // Element type of the dynamic matrix
        , bool SO >   // Storage order of the dynamic matrix
inline Creator< blaze::UniUpperMatrix< blaze::DynamicMatrix<T,SO> > >::Creator( size_t n, const Creator<T>& elementCreator )
   : n_( n )                // The number of rows and columns of the uniupper dynamic matrix
   , ec_( elementCreator )  // Creator for the elements of the uniupper dynamic matrix
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created uniupper dynamic matrix.
//
// \return The randomly generated uniupper dynamic matrix.
*/
template< typename T  // Element type of the dynamic matrix
        , bool SO >   // Storage order of the dynamic matrix
inline blaze::UniUpperMatrix< blaze::DynamicMatrix<T,SO> >
   Creator< blaze::UniUpperMatrix< blaze::DynamicMatrix<T,SO> > >::operator()() const
{
   return (*this)( Default() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a randomly created uniupper dynamic matrix.
//
// \param policy The creation policy for the elements of fundamental data type.
// \return The randomly generated uniupper dynamic matrix.
*/
template< typename T     // Element type of the dynamic matrix
        , bool SO >      // Storage order of the dynamic matrix
template< typename CP >  // Creation policy
inline blaze::UniUpperMatrix< blaze::DynamicMatrix<T,SO> >
   Creator< blaze::UniUpperMatrix< blaze::DynamicMatrix<T,SO> > >::operator()( const CP& policy ) const
{
   blaze::UniUpperMatrix< blaze::DynamicMatrix<T,SO> > matrix( n_ );

   // Initialization of a column-major matrix
   if( SO ) {
      for( size_t j=1UL; j<n_; ++j )
         for( size_t i=0UL; i<j; ++i )
            matrix(i,j) = ec_( policy );
   }

   // Initialization of a row-major matrix
   else {
      for( size_t i=0UL; i<n_; ++i )
         for( size_t j=i+1UL; j<n_; ++j )
            matrix(i,j) = ec_( policy );
   }

   return matrix;
}
//*************************************************************************************************

} // namespace blazetest

#endif
