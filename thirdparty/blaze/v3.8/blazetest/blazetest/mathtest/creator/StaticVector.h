//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/StaticVector.h
//  \brief Specialization of the Creator class template for StaticVector
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_STATICVECTOR_H_
#define _BLAZETEST_MATHTEST_CREATOR_STATICVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/StaticVector.h>
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
/*!\brief Specialization of the Creator class template for static vectors.
//
// This specialization of the Creator class template is able to create random static vectors.
*/
template< typename T  // Element type of the static vector
        , size_t N    // Number of elements of the static vector
        , bool TF >   // Transpose flag of the static vector
class Creator< blaze::StaticVector<T,N,TF> >
{
 public:
   //**Type definitions****************************************************************************
   using Type = blaze::StaticVector<T,N,TF>;  //!< Type to be created by the Creator.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Creator( const Creator<T>& elementCreator = Creator<T>() );
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

   blaze::StaticVector<T,N,TF> operator()() const;

   template< typename CP >
   blaze::StaticVector<T,N,TF> operator()( const CP& policy ) const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Creator<T> ec_;  //!< Creator for the elements of the static vector.
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
/*!\brief Constructor for the creator specialization for StaticVector.
//
// \param elementCreator The creator for the elements of the static vector.
*/
template< typename T  // Element type of the static vector
        , size_t N    // Number of elements of the static vector
        , bool TF >   // Transpose flag of the static vector
inline Creator< blaze::StaticVector<T,N,TF> >::Creator( const Creator<T>& elementCreator )
   : ec_( elementCreator )  // Creator for the elements of the static vector
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created static vector.
//
// \return The randomly generated static vector.
*/
template< typename T  // Element type of the static vector
        , size_t N    // Number of elements of the static vector
        , bool TF >   // Transpose flag of the static vector
inline blaze::StaticVector<T,N,TF> Creator< blaze::StaticVector<T,N,TF> >::operator()() const
{
   return (*this)( Default() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a randomly created static vector.
//
// \param policy The creation policy for the elements of fundamental data type.
// \return The randomly generated static vector.
*/
template< typename T     // Element type of the static vector
        , size_t N       // Number of elements of the static vector
        , bool TF >      // Transpose flag of the static vector
template< typename CP >  // Creation policy
inline blaze::StaticVector<T,N,TF>
   Creator< blaze::StaticVector<T,N,TF> >::operator()( const CP& policy ) const
{
   blaze::StaticVector<T,N,TF> vector;
   for( size_t i=0UL; i<N; ++i )
      vector[i] = ec_( policy );
   return vector;
}
//*************************************************************************************************

} // namespace blazetest

#endif
