//=================================================================================================
/*!
//  \file blazemark/util/MatrixStructure.h
//  \brief Header file for the matrix structure flags
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

#ifndef _BLAZEMARK_UTIL_MATRIXSTRUCTURE_H_
#define _BLAZEMARK_UTIL_MATRIXSTRUCTURE_H_


namespace blazemark {

//=================================================================================================
//
//  MATRIX STRUCTURE FLAGS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Matrix structure flags.
//
// The MatrixStructure enumeration represents all possible structures of (sparse) matrices.
*/
enum MatrixStructure {
   band   = 1,  //!< Flag for banded matrices.
                /*!< The \a band flag indicates a banded matrix of a certain bandwidth. The
                     following example demonstrates a typical 3-banded matrix resulting from
                     a physical setting:
                     \f[\left(\begin{array}{*{5}{c}}
                     2  & -1 & 0  & 0  & 0  \\
                     -1 & 2  & -1 & 0  & 0  \\
                     0  & -1 & 2  & -1 & 0  \\
                     0  & 0  & -1 & 2  & -1 \\
                     0  & 0  & 0  & -1 & 2  \\
                     \end{array}\right)\f]. */
   random = 2   //!< Flag for random matrices.
                /*!< The \a random flag indicates a matrix with randomly determined non-zero
                     entries. The following example demonstrates a random matrix with 2
                     non-zero entries per row:
                     \f[\left(\begin{array}{*{5}{c}}
                     0 & 3 & 0 & 0 & 2 \\
                     1 & 0 & 5 & 0 & 0 \\
                     0 & 1 & 7 & 0 & 0 \\
                     0 & 0 & 1 & 0 & 4 \\
                     0 & 8 & 0 & 1 & 0 \\
                     \end{array}\right)\f]. */
};
//*************************************************************************************************

} // namespace blazemark

#endif
