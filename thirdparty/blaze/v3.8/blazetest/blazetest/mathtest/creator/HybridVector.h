//=================================================================================================
/*!
//  \file blazetest/mathtest/creator/HybridVector.h
//  \brief Specialization of the Creator class template for HybridVector
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

#ifndef _BLAZETEST_MATHTEST_CREATOR_HYBRIDVECTOR_H_
#define _BLAZETEST_MATHTEST_CREATOR_HYBRIDVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/HybridVector.h>
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
/*!\brief Specialization of the Creator class template for hybrid vectors.
//
// This specialization of the Creator class template is able to create random hybrid vectors.
*/
template< typename T  // Element type of the hybrid vector
        , size_t N    // Number of elements of the hybrid vector
        , bool TF >   // Transpose flag of the hybrid vector
class Creator< blaze::HybridVector<T,N,TF> >
{
 public:
   //**Type definitions****************************************************************************
   using Type = blaze::HybridVector<T,N,TF>;  //!< Type to be created by the Creator.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline Creator( const Creator<T>& elementCreator = Creator<T>() );
   explicit inline Creator( size_t size, const Creator<T>& elementCreator = Creator<T>() );
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

   blaze::HybridVector<T,N,TF> operator()() const;

   template< typename CP >
   blaze::HybridVector<T,N,TF> operator()( const CP& policy ) const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t size_;    //!< The size for the hybrid vector.
   Creator<T> ec_;  //!< Creator for the elements of the hybrid vector.
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
/*!\brief Constructor for the creator specialization for HybridVector.
//
// \param elementCreator The creator for the elements of the hybrid vector.
*/
template< typename T  // Element type of the hybrid vector
        , size_t N    // Number of elements of the hybrid vector
        , bool TF >   // Transpose flag of the hybrid vector
inline Creator< blaze::HybridVector<T,N,TF> >::Creator( const Creator<T>& elementCreator )
   : size_( N )             // The size for the hybrid vector
   , ec_( elementCreator )  // Creator for the elements of the hybrid vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the creator specialization for HybridVector.
//
// \param size The size for the hybrid vector.
// \param elementCreator The creator for the elements of the hybrid vector.
*/
template< typename T  // Element type of the hybrid vector
        , size_t N    // Number of elements of the hybrid vector
        , bool TF >   // Transpose flag of the hybrid vector
inline Creator< blaze::HybridVector<T,N,TF> >::Creator( size_t size, const Creator<T>& elementCreator )
   : size_( size )          // The size for the hybrid vector
   , ec_( elementCreator )  // Creator for the elements of the hybrid vector
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns a randomly created hybrid vector.
//
// \return The randomly generated hybrid vector.
*/
template< typename T  // Element type of the hybrid vector
        , size_t N    // Number of elements of the hybrid vector
        , bool TF >   // Transpose flag of the hybrid vector
inline blaze::HybridVector<T,N,TF> Creator< blaze::HybridVector<T,N,TF> >::operator()() const
{
   return (*this)( Default() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a randomly created hybrid vector.
//
// \param policy The creation policy for the elements of fundamental data type.
// \return The randomly generated hybrid vector.
*/
template< typename T     // Element type of the hybrid vector
        , size_t N       // Number of elements of the hybrid vector
        , bool TF >      // Transpose flag of the hybrid vector
template< typename CP >  // Creation policy
inline blaze::HybridVector<T,N,TF>
   Creator< blaze::HybridVector<T,N,TF> >::operator()( const CP& policy ) const
{
   blaze::HybridVector<T,N,TF> vector( size_ );
   for( size_t i=0UL; i<size_; ++i )
      vector[i] = ec_( policy );
   return vector;
}
//*************************************************************************************************

} // namespace blazetest

#endif
