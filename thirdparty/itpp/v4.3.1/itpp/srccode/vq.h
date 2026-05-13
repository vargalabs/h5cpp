/*!
 * \file
 * \brief Definition of a vector quantizer class (unconstrained)
 * \author Thomas Eriksson
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#ifndef VQ_H
#define VQ_H

#include <itpp/base/vec.h>
#include <itpp/base/sort.h>
#include <itpp/base/math/log_exp.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \addtogroup sourcecoding
*/


/*!
  \ingroup sourcecoding
  \brief Class for vector quantization

  The following code illustrates how the VQ can be initialized from a file and
  used to quantize a random vector.
  \code
  Vector_Quantizer Quantizer;
  vec     x,y;
  int     i;

  Quantizer.load("randomvq.vq");
  x=randn(Quantizer.dim());
  i=Quantizer.encode(x);
  y=Quantizer.decode(i);
  \endcode
*/
class ITPP_EXPORT Vector_Quantizer
{
public:
  //! Default constructor
  Vector_Quantizer();
  //! Create a VQ from a VQ file
  Vector_Quantizer(const char *Name);
  //! Encode the input vector
  int encode(const vec &x);
  //! Encode the input vector, and return the num best indices
  ivec encode(const vec &x, int num);
  //! Decode the index
  vec decode(int Index) const;
  //! Decode the indices
  Array<vec> decode(const ivec &Index) const;
  //! Quantize the input vector
  vec Q(const vec &x);
  //! Quantize the input vector
  vec operator()(const vec &x);
  //! Initialize the codebook by a matrix
  void set_codebook(const mat &CB);
  //! Returns the codebook
  mat get_codebook() const;
  //! Set a codevector in the codebook
  void set_codevector(int Index, const vec &indata);
  //! Returns the codevector at the given index
  vec get_codevector(int Index) const;
  //! Rescale and translate a codevector
  void modify_codevector(int no, double mul, const vec &add);
  //! Returns the size (number of codevectors) of the VQ
  int size() const;
  //! Returns the dimension of the VQ
  int dim() const;
  //! Returns the number of bits of the VQ [log2(size)/dim]
  int nobits() const;
  /*!
    \brief Load the codebook from a file
    \param Name The name of the VQ file

    The file format is a textfile where each row is a vector from the codebook.
  */
  void load(const char *Name);
  /*!
    \brief Save the codebook to a file
    \param Name The name of the VQ file

    The file format is a textfile where each row is a vector from the codebook.
  */
  void save(const char *Name) const;
  //! Returns the distortion at the latest time a vector was encoded
  double latest_distortion();
protected:
  //! The vector containing the code book
  vec CodeBook;
  //! The size of the code book
  int Size;
  //! The dimension of the code book
  int Dim;
  //! The distortion at the latest time a vector was encoded
  double LatestDist;
};

// INLINE FUNCTIONS

inline int Vector_Quantizer::size() const { return Size; }
inline int Vector_Quantizer::nobits() const { return levels2bits(Size); }
inline int Vector_Quantizer::dim() const { return Dim; }
inline double Vector_Quantizer::latest_distortion() { return LatestDist; }
inline vec Vector_Quantizer::decode(int Index) const { return get_codevector(Index); }
inline vec Vector_Quantizer::Q(const vec &x) { return decode(encode(x)); }
inline vec Vector_Quantizer::operator()(const vec &x) { return Q(x); }

/*!
  \ingroup sourcecoding
  \brief Class for vector quantization

  The following code illustrates how the quantizer can be initialized from a file and
  used to quantize a random vector.
  \code
  Scalar_Quantizer Quantizer;
  double    x,y;
  int     i;

  Quantizer.load("random.sq");
  x=randn();
  i=Quantizer.encode(x);
  y=Quantizer.decode(i);
  \endcode
*/

class ITPP_EXPORT Scalar_Quantizer
{
public:
  //! Default constructor
  Scalar_Quantizer();
  //! Create a VQ from a VQ file
  Scalar_Quantizer(const char *Name);
  //! Encode
  int encode(double x) const;
  //! Encode the input vector
  ivec encode(const vec &x) const;
  //! Decode the index
  double decode(int Index) const;
  //! Decode the indices
  vec decode(const ivec &Index) const;
  //! Quantize
  double Q(double x) const;
  //! Quantize the input vector
  vec Q(const vec &x) const;
  //! Quantize
  double operator()(double x) const;
  //! Quantize the input vector
  vec operator()(const vec &x) const;
  //! Initialize the codebook by a matrix
  void set_levels(const vec &L);
  //! Returns the codebook
  vec get_levels() const;
  //! Returns the size (number of codevectors) of the VQ
  int size() const;
protected:
  //! The vector containing the code book
  vec Levels;
  //! The distortion at the latest time a vector was encoded
  double LatestDist;
};

inline int Scalar_Quantizer::size() const { return Levels.length(); }
inline double Scalar_Quantizer::decode(int Index) const { return Levels(Index); }
inline double Scalar_Quantizer::Q(double x) const { return decode(encode(x)); }
inline double Scalar_Quantizer::operator()(double x) const { return Q(x); }
inline vec Scalar_Quantizer::operator()(const vec &x)  const { return Q(x); }
inline void Scalar_Quantizer::set_levels(const vec &L) {Levels = L;sort(Levels); }
inline vec Scalar_Quantizer::get_levels() const {return Levels; }

//! ADD DOCUMENTATION HERE
ITPP_EXPORT int scalar_encode(double x, vec &Levels) ;
//! ADD DOCUMENTATION HERE
ITPP_EXPORT ivec scalar_encode(vec &x, vec &Levels);
//! ADD DOCUMENTATION HERE
inline double scalar_quantize(double x, vec &Levels) { return Levels(scalar_encode(x, Levels)); }
//! ADD DOCUMENTATION HERE
inline vec scalar_quantize(vec &x, vec &Levels) { return Levels(scalar_encode(x, Levels)); }


} // namespace itpp

#endif // #ifndef VQ_H
