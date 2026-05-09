/*!
 * \file
 * \brief Definitions of a Hamming code class
 * \author Tony Ottosson
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

#ifndef HAMMING_H
#define HAMMING_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/comm/channel_code.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \ingroup fec
  \brief Binary Hamming codes
*/
class ITPP_EXPORT Hamming_Code : public Channel_Code
{
public:
  //! Constructor for \c hamming(n,k). n = pow(2,m)-1 and k = pow(2,m)-m-1.
  Hamming_Code(int m);

  //! Destructor
  virtual ~Hamming_Code() { }

  //! Hamming encoder. Will truncate some bits if not \a length = \c integer * \a k.
  virtual void encode(const bvec &uncoded_bits, bvec &coded_bits);
  //! Hamming encoder. Will truncate some bits if not \a length = \c integer * \a k.
  virtual bvec encode(const bvec &uncoded_bits);

  //! Hamming decoder. Will truncate some bits if not \a length = \c integer * \a n.
  virtual void decode(const bvec &coded_bits, bvec &decoded_bits);
  //! Hamming decoder. Will truncate some bits if not \a length = \c integer * \a n.
  virtual bvec decode(const bvec &coded_bits);

  // Soft-decision decoding is not implemented
  virtual void decode(const vec &received_signal, bvec &output);
  virtual bvec decode(const vec &received_signal);

  //! Get the code rate
  virtual double get_rate() const { return static_cast<double>(k) / n; };

  //! Gets the code length \a n.
  int get_n() const { return n; };
  //! Gets the number of information bits per code word, \a k.
  int get_k() const { return k; };
  //! Gets the parity check matrix for the code.
  bmat get_H() const { return H; };
  //! Gets the generator matrix for the code.
  bmat get_G() const { return G; };
private:
  int n, k;
  bmat H, G;
  void generate_H(void);
  void generate_G(void);
};

} // namespace itpp

#endif // #ifndef HAMMING_H
