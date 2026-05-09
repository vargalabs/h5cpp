/*!
 * \file
 * \brief Definition of the Extended Golay Code (24, 12, 8)
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

#ifndef EGOLAY_H
#define EGOLAY_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/comm/channel_code.h>
#include <itpp/itexports.h>

namespace itpp
{

/*!
  \ingroup fec
  \brief Extended Golay code (24,12,8).
  \author Tony Ottosson

  The code is given in systematic form with the information bits
  first, followed by the parity check bits. The decoder uses the
  arithmetic decoding algorithm that is for example described in
  Wicker "Error Control Systems for Digital Communication and
  Storage", Prentice Hall, 1995 (page 143).
*/
class ITPP_EXPORT Extended_Golay : public Channel_Code
{
public:
  //! Constructor
  Extended_Golay();
  //! Destructor
  virtual ~Extended_Golay() { }

  //! Encoder. Will truncate some bits if not \a length = \c integer * 12
  virtual void encode(const bvec &uncoded_bits, bvec &coded_bits);
  //! Encoder. Will truncate some bits if not \a length = \c integer * 12
  virtual bvec encode(const bvec &uncoded_bits);

  //! Decoder. Will truncate some bits if not \a length = \c integer * 24
  virtual void decode(const bvec &coded_bits, bvec &decoded_bits);
  //! Decoder. Will truncate some bits if not \a length = \c integer * 24
  virtual bvec decode(const bvec &coded_bits);

  // Soft-decision decoding is not implemented
  virtual void decode(const vec &received_signal, bvec &output);
  virtual bvec decode(const vec &received_signal);

  //! Get the code rate
  virtual double get_rate() const { return 0.5; };

  //! Gets the generator matrix for the code (also the parity check matrix)
  bmat get_G() const { return G; }
private:
  bmat B, G;
};

} // namespace itpp

#endif // #ifndef EGOLAY_H
